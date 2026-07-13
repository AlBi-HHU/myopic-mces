// ── One-dataset MCES workflow ────────────────────────────────────────────────
// All pairwise MCES distances within a single SMILES list.
//
// Usage:
//   nextflow run nextflow/one_dataset.nf --smiles smiles.txt
//
// Options:
//   --threshold 10 --batch_size 224960000 --use_db_lookup --out results

include { PREPARE_INPUT_ONE; COMPUTE_BATCH; COMBINE_BATCHES; SANITY_CHECK } from './modules'

workflow {
    if (!params.smiles) {
        exit 1, "Provide --smiles <file.txt> (one SMILES per line)"
    }

    // 1. Prepare
    batches_ch = PREPARE_INPUT_ONE(params.smiles).flatten()

    // Split: precomputed batch (from db lookup) vs batches to compute
    precomputed_ch = batches_ch.filter { it.name ==~ /batch0_precomputed\.hdf5/ }
    compute_ch    = batches_ch.filter { it.name ==~ /batch\d+\.hdf5/ }

    // 2. Compute (parallel, one task per batch)
    computed_ch = COMPUTE_BATCH(compute_ch)

    // 3. Combine (computed + precomputed)
    all_batches = computed_ch.mix(precomputed_ch).collect()
    combined_ch = COMBINE_BATCHES(all_batches, '')

    // 4. Sanity check
    if (params.sanity_check) {
        SANITY_CHECK(combined_ch, file('nextflow/sanity_check.py'))
    }
}
