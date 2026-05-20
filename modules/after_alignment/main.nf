process after_alignment {
    publishDir "${params.outDir}/results", mode: 'copy'

    input:
    path inputDir

    output:
    path "concatenated__allgenes_allspecies.fasta"
    path "genes_length_after.tsv"
    path "partition.nex"
    path "partition.raxml"

    script:
    """
    species_list=( ${params.species.collect { '"' + it + '"' }.join(' ')} )
    mito_genes_list=( ${params.mito_genes.collect { '"' + it + '"' }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { '"' + it + '"' }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )

    output_file="concatenated__allgenes_allspecies.fasta"
    > "\$output_file"

    # Collect gene lengths first
    declare -A gene_lengths
    for gene in "\${all_genes_list[@]}"; do
        length=\$(awk '
        /^>/ { if (seq) exit; next }
        {
            gsub(/[[:space:]]/, "", \$0)
            seq += length(\$0)
        }
        END { print seq }
        ' "${inputDir}/concatenated_\${gene}_aligned.fas")
        gene_lengths[\$gene]=\$length
    done

    for species in "\${species_list[@]}"; do
        species_nospace=\$(echo "\$species" | tr ' ' '_')
        concatenated_sequence=""
        for gene in "\${all_genes_list[@]}"; do
            aligned_file="${inputDir}/concatenated_\${gene}_aligned.fas"
            gene_sequence=\$(awk -v species="\$species" '
            BEGIN {found=0}
            /^>/ {
                found = (\$0 ~ ("organism=" species))
                next
            }
            found {
                gsub(/[[:space:]]/, "", \$0)
                printf "%s", \$0
            }
            ' "\$aligned_file")
            
            # If gene sequence not found, add gap characters
            if [ -z "\$gene_sequence" ]; then
                gene_len=\${gene_lengths[\$gene]:-0}
                if [ "\$gene_len" -gt 0 ]; then
                    gene_sequence=\$(printf -- '-%.0s' \$(seq 1 \$gene_len))
                fi
            fi
            concatenated_sequence="\${concatenated_sequence}\${gene_sequence}"
        done
        echo ">\$species_nospace" >> "\$output_file"
        printf '%s\n' "\$concatenated_sequence" >> "\$output_file"
    done

    > genes_length_after.tsv
    for gene in "\${all_genes_list[@]}"; do
        # Count aligned columns, keeping gap characters in the partition span.
        length=\$(awk '
        /^>/ { if (seq) exit; next }
        {
            gsub(/[[:space:]]/, "", \$0)
            seq += length(\$0)
        }
        END { print seq }
        ' "${inputDir}/concatenated_\${gene}_aligned.fas")
        echo -e "\${gene}\t\${length}" >> genes_length_after.tsv
    done

    # Generate Nexus partition file
    partition_file="partition.nex"
    {
        echo "#nexus"
        echo "begin sets;"
        
        # Add mitochondrial genes as charsets
        echo "    \"Mitochondrial genes:\" "
        start=1
        for gene in "\${mito_genes_list[@]}"; do
            length=\$(grep "^\${gene}\t" genes_length_after.tsv | cut -f2)
            if [ -n "\$length" ]; then
                end=\$((start + length - 1))
                echo "    charset \${gene} = \${start}-\${end};"
                start=\$((end + 1))
            fi
        done
        
        # Add nuclear genes as charsets
        echo "    \"Nuclear genes:\" "
        for gene in "\${nuclear_genes_list[@]}"; do
            length=\$(grep "^\${gene}\t" genes_length_after.tsv | cut -f2)
            if [ -n "\$length" ]; then
                end=\$((start + length - 1))
                echo "    charset \${gene} = \${start}-\${end};"
                start=\$((end + 1))
            fi
        done
        
        echo "end;"
    } > "\${partition_file}"

    # Generate a RAxML-style partition file for IQ-TREE compatibility
    partition_raxml_file="partition.raxml"
    > "\${partition_raxml_file}"
    start=1
    for gene in "\${all_genes_list[@]}"; do
        length=\$(grep "^\${gene}\t" genes_length_after.tsv | cut -f2)
        if [ -n "\$length" ]; then
            end=\$((start + length - 1))
            echo "DNA, \${gene} = \${start}-\${end}" >> "\${partition_raxml_file}"
            start=\$((end + 1))
        fi
    done
    """
}
