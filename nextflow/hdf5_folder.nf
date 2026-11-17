// ── HDF5-folder MCES workflow ────────────────────────────────────────────────
// A folder of many small one-dataset HDF5 files (SMILES in 'smiles' or
// 'mces_smiles_order'). Files without 'mces' are aggregated into compute
// batches; results are written back into the input files in place
// ('mces', 'computation_modes_order', 'computation_args').
// With --sanity_check (default), one randomly picked written file is re-checked.
//
// Usage:
//   nextflow run nextflow/hdf5_folder.nf --hdf5_dir /path/to/folder

include { PACK_HDF5_FOLDER; COMPUTE_BATCH; UNPACK_HDF5_FOLDER; SANITY_CHECK } from './modules'

workflow HDF5_FOLDER {
    if (!params.hdf5_dir) {
        exit 1, "Provide --hdf5_dir <folder with .hdf5 files>"
    }
    batches_ch = PACK_HDF5_FOLDER(file(params.hdf5_dir).toAbsolutePath().toString()).flatten()
    written_ch = UNPACK_HDF5_FOLDER(COMPUTE_BATCH(batches_ch))

    // Sanity check on one randomly picked file written in this run
    if (params.sanity_check) {
        sample_ch = written_ch
            .flatMap { it.readLines().findAll { l -> l.trim() } }
            .toList()
            .filter { it }
            .map { files -> file(files[new Random().nextInt(files.size())]) }
            .view { "sanity check on: ${it}" }
        SANITY_CHECK(sample_ch, file("${moduleDir}/sanity_check.py"))
    }
}

workflow {
    HDF5_FOLDER()
}
