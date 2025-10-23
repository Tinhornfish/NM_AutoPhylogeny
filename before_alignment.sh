#!/usr/bin/env bash

mkdir -p data

species_list=(
    "Pan paniscus"
    "Pan troglodytes"
    "Gorilla beringei"
    "Gorilla gorilla"
    "Pongo abelii"
    "Pongo pygmaeus"
    "Hoolock hoolock"
    "Hoolock leuconedys"
    "Hylobates agilis"
    "Nomascus leucogenys"
    "Cercocebus atys"
    "Macaca mulatta"
    "Macaca fascicularis"
)
genes_list=(ND4 ND6 COX1 ND2 CYTB COX2 ND5 ATP8 ND3 ND4L ND1 COX3)

touch "data/genes_lengths_before.tsv"

for species in "${species_list[@]}"; do
    species_nospace=$(echo "$species" | tr ' ' '_')
    mkdir -p "data/$species_nospace"
    for gene in "${genes_list[@]}"; do
        echo "Downloading gene $gene for species $species..."
        datasets download gene symbol "$gene" \
            --taxon "$species" \
            --include gene \
            --filename "data/$species_nospace/${gene}.zip"
        
        if [ $? -ne 0 ]; then
            echo "Error downloading gene $gene for species $species. Skipping..."
            continue
        fi

        unzip "data/$species_nospace/${gene}.zip" -d "data/$species_nospace/${gene}"
        if [ $? -ne 0 ]; then
            echo "Error unzipping data for gene $gene for species $species. Skipping..."
            continue
        fi

        gene_length=$(grep -v ">" "data/$species_nospace/${gene}/ncbi_dataset/data/gene.fna" | wc -c)
        if [ $? -ne 0 ]; then
            echo "Error calculating gene length for $gene for species $species. Skipping..."
            continue
        fi

        echo -e "$species_nospace\t$gene\t$gene_length" >> "data/genes_lengths_before.tsv"
    done
done

for gene in "${genes_list[@]}"; do
    touch "data/concatenated_${gene}.fasta"
    for species in "${species_list[@]}"; do
        species_nospace=$(echo "$species" | tr ' ' '_')
        cat "data/$species_nospace/${gene}/ncbi_dataset/data/gene.fna" >> "data/concatenated_${gene}.fasta"
    done
done
