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

output_file="concatenated__allgenes_allspecies.fasta"
touch "$output_file"  # Create the output file if it doesn't exist
> "$output_file"  # Clear the output file if it exists

for species in "${species_list[@]}"; do
    species_nospace=$(echo "$species" | tr ' ' '_')
    echo ">$species_nospace" >> "$output_file"
    for gene in "${genes_list[@]}"; do
        awk -v species="$species" '
        BEGIN {found=0}
        $0 ~ "organism="species {found=1; next}
        found && $0 ~ /^>/ {found=0}
        found {print}
        ' "data/concatenated_${gene}_aligned.fas" >> "$output_file"
    done
    echo >> "$output_file"  # Add a newline after each species
done

touch genes_length_after.tsv
for gene in "${genes_list[@]}"; do
    length=$(awk '/^>/ {if (seq) exit; next} {seq += length($0)} END {print seq}' "data/concatenated_${gene}_aligned.fas")
    echo -e "${gene}\t${length}" >> genes_length_after.tsv
done
