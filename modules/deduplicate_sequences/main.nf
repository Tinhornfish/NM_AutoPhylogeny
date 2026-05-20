process deduplicate_sequences {
    publishDir "${params.outDir}", mode: 'copy'

    input:
    path inputDir

    output:
    path "data"

    script:
    """
    set -e

    mkdir -p dedup_out
    cp -r "${inputDir}"/. dedup_out/

    mito_genes_list=( ${params.mito_genes.collect { "'" + it + "'" }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { "'" + it + "'" }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )

    for gene in "\${all_genes_list[@]}"; do
        in_file="dedup_out/concatenated_\${gene}.fasta"
        out_file="dedup_out/concatenated_\${gene}.dedup.fasta"

        if [ ! -f "\$in_file" ]; then
            echo "Warning: \$in_file not found, skipping"
            continue
        fi

        python3 "${workflow.projectDir}/bin/dedup.py" "\$in_file" "\$out_file"
        mv "\$out_file" "\$in_file"
    done

    rm -rf data
    mv dedup_out data
    """
}
