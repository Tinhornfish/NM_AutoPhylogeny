process alignment {
    publishDir "${params.outDir}", mode: 'copy'

    container 'thetinhornfish/clustalo:1.2.4'

    input:
    path inputDir

    output:
    path "data"

    script:
    """
    set -e

    # Build bash arrays from Nextflow params
    mito_genes_list=( ${params.mito_genes.collect { "'" + it + "'" }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { "'" + it + "'" }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )

    for gene in "\${all_genes_list[@]}"; do
        in_file="${inputDir}/concatenated_\${gene}.fasta"
        out_file="data/concatenated_\${gene}_aligned.fas"
        if [ -f "\$in_file" ]; then
            echo "Aligning \$in_file -> \$out_file"
            clustalo -i "\$in_file" -o "\$out_file" --force
        else
            echo "Warning: input file \$in_file not found, skipping"
        fi
    done
    """
  
}
