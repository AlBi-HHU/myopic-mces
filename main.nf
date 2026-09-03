// ── Main entry point ─────────────────────────────────────────────────────────
// Dispatches to the one-dataset or two-dataset MCES workflow based on params.
//
// Usage (local checkout):
//   nextflow run . --smiles smiles.txt                       // one-dataset
//   nextflow run . --smiles_a a.txt --smiles_b b.txt          // two-dataset
//
// Usage (directly from GitHub):
//   nextflow run AlBi-HHU/myopic-mces -r <branch> --smiles smiles.txt
//
// Options: --threshold 10 --batch_size 224960000 --use_db_lookup --out results

include { ONE_DATASET } from './nextflow/one_dataset'
include { TWO_DATASETS } from './nextflow/two_datasets'

workflow {
    if (params.smiles_a || params.smiles_b) {
        if (!params.smiles_a || !params.smiles_b) {
            exit 1, "Two-dataset mode needs both --smiles_a <file.txt> and --smiles_b <file.txt>"
        }
        TWO_DATASETS()
    } else if (params.smiles) {
        ONE_DATASET()
    } else {
        exit 1, "Provide --smiles <file.txt> (one-dataset) or --smiles_a/--smiles_b <file.txt> (two-dataset)"
    }
}
