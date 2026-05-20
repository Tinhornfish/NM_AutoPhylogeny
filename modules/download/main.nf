process download {
    publishDir "${params.outDir}", mode: 'copy'

    container 'deptmetagenom/ncbi-datasets-cli:18.9.0'

    input:
    val species
    val mito_genes
    val nuclear_genes

    output:
    path "data"

    script:
    """
    set -e

    mkdir -p data

    species_list=( ${species.collect { '"' + it + '"' }.join(' ')} )
    mito_genes_list=( ${mito_genes.collect { '"' + it + '"' }.join(' ')} )
    nuclear_genes_list=( ${nuclear_genes.collect { '"' + it + '"' }.join(' ')} )

    touch "data/genes_lengths_before.tsv"
    touch "data/missing_queries.tsv"

    for species in "\${species_list[@]}"; do
        species_nospace=\$(echo "\$species" | tr ' ' '_')
        mkdir -p "data/\$species_nospace"

        for gene in "\${mito_genes_list[@]}"; do
            echo "Downloading mito gene \$gene for \$species..."
            if ! datasets download gene symbol "\$gene" \
                --taxon "\$species" \
                --include gene \
                --filename "data/\$species_nospace/\${gene}.zip"; then
                echo "datasets download failed for \$gene / \$species; trying eutils fallback anyway..."
            fi
            if [ -s "data/\$species_nospace/\${gene}.zip" ]; then
                if ! unzip -o "data/\$species_nospace/\${gene}.zip" -d "data/\$species_nospace/\${gene}"; then
                    echo "unzip failed for \$gene / \$species"
                fi
            fi
            gene_file="data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna"
            gene_dir="data/\$species_nospace/\${gene}/ncbi_dataset/data"
            if [ -s "\$gene_file" ]; then
                tmp_file="\${gene_file}.tmp"
                awk -v species="\$species" -v gene="\$gene" '
                    /^>/ { print \$0 " [organism=" species "] [gene=" gene "]"; next }
                    { print }
                ' "\$gene_file" > "\$tmp_file"
                mv "\$tmp_file" "\$gene_file"

                gene_length=\$(awk '
                    /^>/ { next }
                    { gsub(/[[:space:]]/, "", \$0); n += length(\$0) }
                    END { print n }
                ' "\$gene_file")

                printf '%s\t%s\t%s\n' "\$species_nospace" "\$gene" "\$gene_length" >> "data/genes_lengths_before.tsv"
            else
                mkdir -p "\$gene_dir"
                printf '%s\t%s\n' "\$species_nospace" "\$gene" >> "data/missing_queries.tsv"
            fi
        done
            for gene in "\${nuclear_genes_list[@]}"; do
                echo "Downloading nuclear gene \$gene for \$species..."
                if ! datasets download gene symbol "\$gene" \
                    --taxon "\$species" \
                    --include gene \
                    --filename "data/\$species_nospace/\${gene}.zip"; then
                    echo "datasets download failed for \$gene / \$species; trying eutils fallback anyway..."
                fi
                if [ -s "data/\$species_nospace/\${gene}.zip" ]; then
                    if ! unzip -o "data/\$species_nospace/\${gene}.zip" -d "data/\$species_nospace/\${gene}"; then
                        echo "unzip failed for \$gene / \$species"
                    fi
                fi
                gene_file="data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna"
                gene_dir="data/\$species_nospace/\${gene}/ncbi_dataset/data"
                if [ -s "\$gene_file" ]; then
                    tmp_file="\${gene_file}.tmp"
                    awk -v species="\$species" -v gene="\$gene" '
                        /^>/ { print \$0 " [organism=" species "] [gene=" gene "]"; next }
                        { print }
                    ' "\$gene_file" > "\$tmp_file"
                    mv "\$tmp_file" "\$gene_file"

                    gene_length=\$(awk '
                        /^>/ { next }
                        { gsub(/[[:space:]]/, "", \$0); n += length(\$0) }
                        END { print n }
                    ' "\$gene_file")

                    printf '%s\t%s\t%s\n' "\$species_nospace" "\$gene" "\$gene_length" >> "data/genes_lengths_before.tsv"
                else
                    mkdir -p "\$gene_dir"
                    printf '%s\t%s\n' "\$species_nospace" "\$gene" >> "data/missing_queries.tsv"
                fi
            done
    done

    """
}
