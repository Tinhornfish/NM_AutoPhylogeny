nextflow.enable.dsl=2

params.mito_genes = [
    'CYTB'
]
params.nuclear_genes = [
    'RAG1'
]
params.species = [
    'Canis lupus familiaris',
    'Canis adustus',
    'Canis aureus',
    'Canis latrans',
    'Chrysocyon brachyurus',
    'Lycaon pictus',
    'Cerdocyon thous',
    'Cuon alpinus',
    'Urocyon cinereoargenteus',
    'Nyctereutes procyonoides',
    'Otocyon megalotis',
    'Vulpes lagopus',
    'Vulpes vulpes',
    'Speothos venaticus',
    'Canis simensis',
    'Felis catus'
]
params.outDir = 'Results'

include { download } from './modules/download/main.nf'
include { download_fallback } from './modules/download_fallback/main.nf'
include { concatenate_sequences } from './modules/concatenate_sequences/main.nf'
include { nuccore_efetch_fallback } from './modules/nuccore_efetch_fallback/main.nf'
include { deduplicate_sequences } from './modules/deduplicate_sequences/main.nf'
include { pre_alignment_validation } from './modules/pre_alignment_validation/main.nf'
include { alignment } from './modules/alignment/main.nf'
include { after_alignment } from './modules/after_alignment/main.nf'
include { iqtree } from './modules/iqtree/main.nf'
include { tree_image_creation } from './modules/tree_image_creation/main.nf'

workflow {
    download(params.species, params.mito_genes, params.nuclear_genes)
    download_fallback(download.out)
    nuccore_efetch_fallback(download_fallback.out)
    concatenate_sequences(nuccore_efetch_fallback.out, params.mito_genes, params.nuclear_genes, params.species)
    deduplicate_sequences(concatenate_sequences.out)
    pre_alignment_validation(deduplicate_sequences.out)
    alignment(pre_alignment_validation.out)
    after_alignment(alignment.out)
    iqtree(after_alignment.out[0], after_alignment.out[3])
    tree_image_creation(iqtree.out[0])
}
