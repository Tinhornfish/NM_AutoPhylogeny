process concatenate_sequences {
    publishDir "${params.outDir}", mode: 'copy'

    input:
    path inputDir
    val mito_genes
    val nuclear_genes
    val species

    output:
    path "data"

    script:
    """
    set -e

    source_dir=\$(readlink -f "${inputDir}" 2>/dev/null || true)
    if [ -z "\$source_dir" ] || [ ! -d "\$source_dir" ]; then
        source_dir="${inputDir}"
    fi

    rm -rf data
    mkdir -p data
    cp -r "\$source_dir"/. data/

    mito_genes_list=( ${mito_genes.collect { '"' + it + '"' }.join(' ')} )
    nuclear_genes_list=( ${nuclear_genes.collect { '"' + it + '"' }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )
    species_list=( ${species.collect { '"' + it + '"' }.join(' ')} )

    for gene in "\${all_genes_list[@]}"; do
        > "data/concatenated_\${gene}.fasta"
        for sp in "\${species_list[@]}"; do
            species_nospace=\$(echo "\$sp" | tr ' ' '_')
            gene_file="data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna"
            if [ -s "\$gene_file" ]; then
                cat "\$gene_file" >> "data/concatenated_\${gene}.fasta"
            fi
        done
    done
    """
}
