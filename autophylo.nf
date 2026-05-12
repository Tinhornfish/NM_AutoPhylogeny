params.genes = [
    'COX1','CYTB','COX2','ATP8','COX3'
]
params.species = [
    'Pan paniscus',
    'Pan troglodytes',
    'Gorilla beringei',
    'Gorilla gorilla',
    'Pongo abelii',
    'Pongo pygmaeus',
    'Hoolock hoolock',
    'Hoolock leuconedys',
    'Hylobates agilis',
    'Nomascus leucogenys',
    'Cercocebus atys',
    'Macaca mulatta',
    'Macaca fascicularis'
]
params.outDir = 'results'


/*
 * This process downloads all given genes of the user input species, and prepares them for alignment
 */
process before_alignment {
    publishDir "${params.outDir}", mode: 'copy'

    container 'deptmetagenom/ncbi-datasets-cli:18.9.0'

    input:
    val species
    val genes

    output:
    path "data"

    script:
    """
    set -e

    mkdir -p data

    species_list=( ${species.collect { '"' + it + '"' }.join(' ')} )
    genes_list=( ${genes.collect { '"' + it + '"' }.join(' ')} )

    touch "data/genes_lengths_before.tsv"

    for species in "\${species_list[@]}"; do
        species_nospace=\$(echo "\$species" | tr ' ' '_')
        mkdir -p "data/\$species_nospace"
        for gene in "\${genes_list[@]}"; do
            echo "Downloading gene \$gene for species \$species..."
            datasets download gene symbol "\$gene" \
                --taxon "\$species" \
                --include gene \
                --filename "data/\$species_nospace/\${gene}.zip"
            
            if [ \$? -ne 0 ]; then
                echo "Error downloading gene \$gene for species \$species. Skipping..."
                continue
            fi

            unzip "data/\$species_nospace/\${gene}.zip" -d "data/\$species_nospace/\${gene}"
            if [ \$? -ne 0 ]; then
                echo "Error unzipping data for gene \$gene for species \$species. Skipping..."
                continue
            fi

            gene_length=\$(grep -v ">" "data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna" | wc -c)
            if [ \$? -ne 0 ]; then
                echo "Error calculating gene length for \$gene for species \$species. Skipping..."
                continue
            fi

            echo -e "\$species_nospace\t\$gene\t\$gene_length" >> "data/genes_lengths_before.tsv"
        done
    done

    for gene in "\${genes_list[@]}"; do
        touch "data/concatenated_\${gene}.fasta"
        for species in "\${species_list[@]}"; do
            species_nospace=\$(echo "\$species" | tr ' ' '_')
            cat "data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna" >> "data/concatenated_\${gene}.fasta"
        done
    done
    """

}

/*
 * This process aligns the sequences of all species for each gene
 * 
 */
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
    genes_list=( ${params.genes.collect { "'" + it + "'" }.join(' ')} )

    for gene in "\${genes_list[@]}"; do
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

/*
 * This process concatenates all aligned genes so that it1s ready for tree construction
 */
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
    genes_list=( ${params.genes.collect { '"' + it + '"' }.join(' ')} )

    output_file="concatenated__allgenes_allspecies.fasta"
    > "\$output_file"

    for species in "\${species_list[@]}"; do
        species_nospace=\$(echo "\$species" | tr ' ' '_')
        concatenated_sequence=""
        for gene in "\${genes_list[@]}"; do
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
            concatenated_sequence="\${concatenated_sequence}\${gene_sequence}"
        done
        echo ">\$species_nospace" >> "\$output_file"
        printf '%s\n' "\$concatenated_sequence" >> "\$output_file"
    done

    > genes_length_after.tsv
    for gene in "\${genes_list[@]}"; do
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
        start=1
        while IFS=\$(printf '\\t') read -r gene length; do
            end=\$((start + length - 1))
            echo "    charset \${gene} = \${start}-\${end};"
            start=\$((end + 1))
        done < genes_length_after.tsv
        echo "end;"
    } > "\${partition_file}"

    # Generate a RAxML-style partition file for IQ-TREE compatibility
    partition_raxml_file="partition.raxml"
    > "\${partition_raxml_file}"
    start=1
    while IFS=\$(printf '\t') read -r gene length; do
        end=\$((start + length - 1))
        echo "DNA, \${gene} = \${start}-\${end}" >> "\${partition_raxml_file}"
        start=\$((end + 1))
    done < genes_length_after.tsv
    """
}

process iqtree {
    publishDir "${params.outDir}/results", mode: 'copy'

    container 'thetinhornfish/iqtree2:2.0.7'

    input:
    path alignment_file
    path partition_file

    output:
    path "*.treefile"
    path "*.log"
    path "*.iqtree"
    path "*.nwk"

    script:
    """
    echo "Running IQ-TREE with alignment: ${alignment_file} and partitions: ${partition_file}"
    
    iqtree2 -s "${alignment_file}" \
        -spp "${partition_file}" \
        -m MFP \
        -bb 1000 \
        -nt AUTO \
        -pre iqtree_output
    
    # Copy treefile to .nwk format
    cp iqtree_output.treefile iqtree_output.nwk
    """
}

workflow {
    before_alignment(species=params.species, genes=params.genes)
    alignment(before_alignment.out)
    after_alignment(alignment.out)
    iqtree(after_alignment.out[0], after_alignment.out[3])
}