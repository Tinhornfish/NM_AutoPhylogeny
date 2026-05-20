process pre_alignment_validation {
    publishDir "${params.outDir}", mode: 'copy'

    input:
    path inputDir

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

    # Build bash arrays from Nextflow params
    mito_genes_list=( ${params.mito_genes.collect { "'" + it + "'" }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { "'" + it + "'" }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )
    species_list=( ${params.species.collect { '"' + it + '"' }.join(' ')} )
    expected_species_count=\${#species_list[@]}

    # Validation: Check that each gene has exactly one copy from all species
    echo "Validating gene sequence counts before alignment..."
    for gene in "\${all_genes_list[@]}"; do
        in_file="data/concatenated_\${gene}.fasta"
        if [ -f "\$in_file" ]; then
            unique_count=\$(awk '
                /^>/ {
                    line = \$0
                    # Extract organism value from [organism=XXX]
                    if (index(line, "[organism=") > 0) {
                        start = index(line, "[organism=") + 10
                        end = index(substr(line, start), "]")
                        if (end > 0) {
                            species = substr(line, start, end - 1)
                            if (!(species in seen)) {
                                count++
                                seen[species] = 1
                            }
                        }
                    }
                }
                END {
                    print (count > 0) ? count : 0
                }
            ' "\$in_file")

            if [ "\$unique_count" -ne "\$expected_species_count" ]; then
                echo "ERROR: Gene \$gene has \$unique_count unique species in FASTA, but expected \$expected_species_count"
                echo "Missing or duplicate species detected. Check concatenated_\${gene}.fasta"
                exit 1
            fi
            echo "✓ Validation OK: Gene \$gene has \$unique_count species (expected \$expected_species_count)"
        else
            echo "ERROR: Input file \$in_file not found"
            echo "=== Headers in concatenated_\${gene}.fasta ==="
            grep "^>" "\$in_file" | head -20
            exit 1
        fi
    done

    echo ""
    echo "All validations passed. Proceeding with alignment..."
    echo ""
    """
}