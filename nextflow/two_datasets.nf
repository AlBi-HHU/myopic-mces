// ── Two-dataset MCES workflow ────────────────────────────────────────────────
// All cross MCES distances between two SMILES lists.
//   smiles_a = main dataset (input_file argument to prepare_input)
//   smiles_b = extra dataset (hdf5_extra_input_file argument)
//
// Usage:
//   nextflow run nextflow/two_datasets.nf --smiles_a enveda.txt --smiles_b splitpool.txt
//
// Options:
//   --threshold 10 --batch_size 224960000 --use_db_lookup --out results

include { PREPARE_INPUT_TWO; COMPUTE_BATCH; COMBINE_BATCHES; SANITY_CHECK } from './modules'

workflow {
    if (!params.smiles_a || !params.smiles_b) {
        exit 1, "Provide --smiles_a <file.txt> and --smiles_b <file.txt>"
    }

    n_a = file(params.smiles_a).readLines().size()
    n_b = file(params.smiles_b).readLines().size()

    // 1. Prepare
    batches_ch = PREPARE_INPUT_TWO(params.smiles_a, params.smiles_b).flatten()

    precomputed_ch = batches_ch.filter { it.name ==~ /batch0_precomputed\.hdf5/ }
    compute_ch    = batches_ch.filter { it.name ==~ /batch\d+\.hdf5/ }

    // 2. Compute (parallel, one task per batch)
    computed_ch = COMPUTE_BATCH(compute_ch)

    // 3. Combine — two_datasets_shape is "n_extra n_main" = "n_b n_a"
    all_batches = computed_ch.mix(precomputed_ch).collect()
    combined_ch = COMBINE_BATCHES(all_batches, "${n_b} ${n_a}")

    // 4. Sanity check
    if (params.sanity_check) {
        SANITY_CHECK(combined_ch, file('nextflow/sanity_check.py'))
    }
}
