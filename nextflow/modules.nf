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

    input:
        path batch

    output:
        path batch

    script:
        def cplex = params.cplex_home ? "PATH=${params.cplex_home}/bin/x86-64_linux/:\$PATH" : ''
        def solver = params.cplex_home ? params.solver : 'COIN_CMD'
        def dyn = params.choose_bound_dynamically ? '--choose_bound_dynamically' : ''
        """
        ${cplex} myopic_mces \\
            --hdf5_mode "${batch}" tmpout \\
            --threshold ${params.threshold} \\
            --solver ${solver} \\
            --num_jobs ${task.cpus} \\
            --solver_onethreaded \\
            --solver_no_msg \\
            ${dyn}
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
