// ── Shared processes for MCES Nextflow workflows ────────────────────────────
// Included by one_dataset.nf and two_datasets.nf

// ── 1. PREPARE_INPUT (one-dataset) ───────────────────────────────────────────
// Splits a single SMILES list into batched HDF5 files (all pairwise comparisons).

process PREPARE_INPUT_ONE {
    tag 'one-dataset'

    input:
        path smiles_file

    output:
        path 'out/batch*.hdf5', emit: batches

    script:
        def db = params.use_db_lookup ? "--use_db_lookup --threshold ${params.threshold}" : ''
        """
        python3 -m myopic_mces.prepare_input \\
            --hdf5_mode \\
            --batch_size ${params.batch_size} \\
            ${db} \\
            "${smiles_file}" \\
            out
        """
}

// ── 1. PREPARE_INPUT (two-dataset) ───────────────────────────────────────────
// Splits two SMILES lists into batched HDF5 files (all cross comparisons).

process PREPARE_INPUT_TWO {
    tag 'two-datasets'

    input:
        path smiles_file_a
        path smiles_file_b

    output:
        path 'out/batch*.hdf5', emit: batches

    script:
        def db = params.use_db_lookup ? "--use_db_lookup --threshold ${params.threshold}" : ''
        """
        python3 -m myopic_mces.prepare_input \\
            --hdf5_mode \\
            --hdf5_extra_input_file "${smiles_file_b}" \\
            --batch_size ${params.batch_size} \\
            ${db} \\
            "${smiles_file_a}" \\
            out
        """
}

// ── 2. COMPUTE_BATCH ─────────────────────────────────────────────────────────
// Runs myopic_mces on a single batch HDF5 file (modifies it in place).

process COMPUTE_BATCH {
    tag "${batch.name}"
    cpus params.cpus
    stageInMode 'copy'   // copy input (myopic_mces modifies the hdf5 in place)

    // `docplex config --upgrade` (below) patches files inside the container's
    // own venv, which is otherwise read-only under Singularity/Apptainer. A
    // writable overlay is needed for that; it must exist before the
    // container starts, and it's per-task so parallel COMPUTE_BATCH tasks
    // (maxForks) don't race on a shared one — Nextflow cleans it up along
    // with the rest of the task's work dir.
    beforeScript params.cplex_home ? { "mkdir -p ${task.workDir}/cplex_overlay" } : ''
    containerOptions params.cplex_home ? { "--bind ${params.cplex_home} --overlay ${task.workDir}/cplex_overlay" } : ''

    input:
        path batch

    output:
        path batch

    script:
        def solver = params.cplex_home ? params.solver : 'COIN_CMD'
        def cplex_setup = params.cplex_home ? """
        export CPLEX_HOME=${params.cplex_home}
        export PATH=\$CPLEX_HOME/bin/x86-64_linux/:\$PATH
        docplex config --upgrade \$CPLEX_HOME
        """.stripIndent() : ''
        def dyn = params.choose_bound_dynamically ? '--choose_bound_dynamically' : ''
        def time_limit = params.solver_time_limit_seconds ? "--solver_time_limit_seconds ${params.solver_time_limit_seconds}" : ''
        """
        ${cplex_setup}
        myopic_mces \\
            --hdf5_mode "${batch}" tmpout \\
            --threshold ${params.threshold} \\
            --solver ${solver} \\
            --num_jobs -1 \\
            --solver_onethreaded \\
            --solver_no_msg \\
            ${dyn} \\
            ${time_limit}
        """
}

// ── 3. COMBINE_BATCHES ───────────────────────────────────────────────────────
// Combines all batch HDF5 results into a single file.

process COMBINE_BATCHES {
    tag 'combine'
    publishDir "${params.out}", mode: 'copy', pattern: 'combined.hdf5'

    input:
        path batches
        val two_datasets_shape

    output:
        path 'combined.hdf5', emit: combined

    script:
        def shape = two_datasets_shape ? "--two_datasets_shape ${two_datasets_shape}" : ''
        """
        python3 -m myopic_mces.combine_hdf5_batches \\
            ${batches} \\
            --out combined.hdf5 \\
            ${shape}
        """
}

// ── 4. SANITY_CHECK ──────────────────────────────────────────────────────────
// Recomputes a random sample of pairs and compares to the combined HDF5.
// Auto-detects one-dataset vs two-dataset mode from the HDF5 keys.

process SANITY_CHECK {
    tag 'sanity'
    publishDir "${params.out}", mode: 'copy', pattern: 'sanity_report.txt'

    input:
        path combined
        path script

    output:
        path 'sanity_report.txt', emit: report

    script:
        def always_stronger = !params.choose_bound_dynamically ? 'True' : 'False'
        """
        python3 "${script}" "${combined}" ${params.threshold} ${always_stronger} ${params.sanity_samples}
        """
}

// ── 5. WRITE_TO_DB ───────────────────────────────────────────────────────────
// Upserts the combined HDF5 results into the mces database via the
// `import-mces` CLI from the mces-database package. Connection parameters
// come from the standard PostgreSQL environment variables.

process WRITE_TO_DB {
    tag 'write_to_db'
    publishDir "${params.out}", mode: 'copy', pattern: 'db_import.log'

    input:
        path combined
        val format

    output:
        path 'db_import.log', emit: log

    script:
        """
        import-mces "${combined}" --format ${format} > db_import.log 2>&1
        """
}
