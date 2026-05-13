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
params.outDir = 'Mayer_Adam'


/*
 * This process downloads all given genes of the user input species, and prepares them for alignment
 */
process before_alignment {
    publishDir "${params.outDir}", mode: 'copy'

    container 'deptmetagenom/ncbi-datasets-cli:18.9.0'

    input:
    val species

    output:
    path "data"

    script:
    """
    set -e

    mkdir -p data

    # curl with retries + exponential backoff and jitter
    curl_retry() {
        max_attempts=3
        attempt=1
        while true; do
            curl -sG "\$@" && return 0
            if [ "\$attempt" -ge "\$max_attempts" ]; then
                return 1
            fi
            backoff=\$((2 ** (attempt - 1)))
            jitter=\$((RANDOM % 3))
            sleep_time=\$((backoff + jitter))
            sleep \$sleep_time
            attempt=\$((attempt + 1))
        done
    }

    # helper to fetch efetch with retry and try multiple rettypes
    efetch_try() {
        local out_file="\$1"
        url_base="\$2"
        seq_id="\$3"
        sleep \$((RANDOM % 2))

        # first try fasta_cds_na (preferred)
        if curl -fsSL -o "\${out_file}.tmp" "\${url_base}&id=\${seq_id}&rettype=fasta_cds_na&retmode=text"; then
            if [ -s "\${out_file}.tmp" ]; then
                mv "\${out_file}.tmp" "\${out_file}" && return 0
            else
                rm -f "\${out_file}.tmp"
            fi
        fi

        sleep \$((RANDOM % 2))
        # fallback to full fasta
        if curl -fsSL -o "\${out_file}.tmp" "\${url_base}&id=\${seq_id}&rettype=fasta&retmode=text"; then
            if [ -s "\${out_file}.tmp" ]; then
                mv "\${out_file}.tmp" "\${out_file}" && return 0
            else
                rm -f "\${out_file}.tmp"
            fi
        fi

        return 1
    }

    species_list=( ${species.collect { '"' + it + '"' }.join(' ')} )
    mito_genes_list=( ${params.mito_genes.collect { '"' + it + '"' }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { '"' + it + '"' }.join(' ')} )

    touch "data/genes_lengths_before.tsv"
    touch "data/missing_queries.tsv"
    touch "data/efetch_recovered.tsv"

    download_and_normalize() {
        local species="\$1"
        local gene="\$2"
        local species_nospace="\$3"
        local gene_file="data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna"
        local out_dir="data/\$species_nospace/\${gene}/ncbi_dataset/data"
        local out_file="\${out_dir}/gene.fna"

        if [ ! -s "\$gene_file" ]; then
            echo "Missing gene.fna for \$gene / \$species. Trying efetch fallback..."
            mkdir -p "\$out_dir"

            query="\${gene}[Gene Name] AND \${species}[Organism]"
                # run esearch with retries
                if ! search_xml=\$(curl_retry --data-urlencode "db=nuccore" --data-urlencode "retmax=20" --data-urlencode "term=\$query" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" 2>/dev/null); then
                    search_xml=""
                fi

            id_list=\$(printf '%s' "\$search_xml" | grep -oE '<Id>[0-9]+</Id>' | sed 's:<Id>::g; s:</Id>::g')

            selected=0
            if [ -n "\$id_list" ]; then
                while IFS= read -r seq_id; do
                    [ -z "\$seq_id" ] && continue
                    tmp_fasta="\${out_file}.tmp"
                    if ! efetch_try "\$tmp_fasta" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore" "\$seq_id"; then
                        rm -f "\$tmp_fasta" || true
                        continue
                    fi

                    # Filter to only the record matching [gene=GENE]
                    awk -v gene="\$gene" '
                        /^>/ { keep = (tolower(\$0) ~ tolower("\\[gene=" gene "\\]")) }
                        keep { print }
                    ' "\$tmp_fasta" > "\${tmp_fasta}.filtered"

                    if [ -s "\${tmp_fasta}.filtered" ]; then
                        mv "\${tmp_fasta}.filtered" "\$tmp_fasta"
                    else
                        rm -f "\$tmp_fasta" "\${tmp_fasta}.filtered"
                        continue
                    fi

                    seq_len=\$(awk '
                        BEGIN { n=0; seen=0 }
                        /^>/ { if (seen==1) exit; seen=1; next }
                        seen { gsub(/[[:space:]]/, "", \$0); n += length(\$0) }
                        END { print n }
                    ' "\$tmp_fasta")

                    if [ -z "\$seq_len" ] || [ "\$seq_len" -lt 200 ]; then
                        rm -f "\$tmp_fasta"
                        continue
                    fi

                    awk -v species="\$species" -v gene="\$gene" '
                        BEGIN { header_seen=0; done=0 }
                        /^>/ {
                            if (header_seen==1) { done=1; exit }
                            header_seen=1
                            print \$0 " [organism=" species "] [gene=" gene "]"
                            next
                        }
                        { if (!done) print }
                    ' "\$tmp_fasta" > "\$out_file"

                    rm -f "\$tmp_fasta"

                    if [ -s "\$out_file" ]; then
                        printf '%s\t%s\t%s\n' "\$species_nospace" "\$gene" "\$seq_id" >> "data/efetch_recovered.tsv"
                        selected=1
                        break
                    fi
                done <<< "\$id_list"
            fi

            # ── Third-level fallback: species-wide GenBank inspection ──
            if [ "\$selected" -eq 0 ]; then
                echo "No suitable gene-specific nuccore hit for \$species / \$gene; trying species-wide GenBank inspection..."

                if ! species_search_xml=\$(curl_retry --data-urlencode "db=nuccore" --data-urlencode "retmax=200" --data-urlencode "term=\${species}[Organism]" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" 2>/dev/null); then
                    species_search_xml=""
                fi

                species_id_list=\$(printf '%s' "\$species_search_xml" | grep -oE '<Id>[0-9]+</Id>' | sed 's:<Id>::g; s:</Id>::g')

                if [ -n "\$species_id_list" ]; then
                    while IFS= read -r seq_id; do
                        [ -z "\$seq_id" ] && continue

                        gb_tmp="\${out_file}.\${seq_id}.gb"
                        if ! curl_retry -o "\$gb_tmp" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\${seq_id}&rettype=gb&retmode=text" 2>/dev/null; then
                            rm -f "\$gb_tmp" || true
                            continue
                        fi

                        if grep -qiE "/gene=\"\${gene}\"|/product=.*\${gene}|\${gene}" "\$gb_tmp"; then
                            if ! efetch_try "\${out_file}.tmp" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore" "\${seq_id}"; then
                                rm -f "\${out_file}.tmp" "\$gb_tmp" || true
                                continue
                            fi

                            if [ -s "\${out_file}.tmp" ]; then
                                awk -v species="\$species" -v gene="\$gene" '
                                    BEGIN { header_seen=0; done=0 }
                                    /^>/ {
                                        if (header_seen==1) { done=1; exit }
                                        header_seen=1
                                        print \$0 " [organism=" species "] [gene=" gene "]"
                                        next
                                    }
                                    { if (!done) print }
                                ' "\${out_file}.tmp" > "\$out_file"

                                rm -f "\${out_file}.tmp" "\$gb_tmp"

                                if [ -s "\$out_file" ]; then
                                    seq_len=\$(awk '
                                        BEGIN { n=0; seen=0 }
                                        /^>/ { if (seen==1) exit; seen=1; next }
                                        seen { gsub(/[[:space:]]/, "", \$0); n += length(\$0) }
                                        END { print n }
                                    ' "\$out_file")

                                    if [ -n "\$seq_len" ] && [ "\$seq_len" -ge 100 ]; then
                                        printf '%s\t%s\t%s\n' "\$species_nospace" "\$gene" "\$seq_id" >> "data/efetch_recovered.tsv"
                                        selected=1
                                        break
                                    else
                                        rm -f "\$out_file"
                                    fi
                                fi
                            fi
                        fi
                        rm -f "\$gb_tmp"
                    done <<< "\$species_id_list"
                fi

                if [ "\$selected" -eq 0 ]; then
                    echo "No suitable nuccore hit for \$species / \$gene after species-wide search"
                    printf '%s\t%s\n' "\$species_nospace" "\$gene" >> "data/missing_queries.tsv"
                    rm -f "\$out_file"
                    return
                fi
            fi
        fi

        # Normalize headers
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
    }

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
                    echo "unzip failed for \$gene / \$species; trying eutils fallback anyway..."
                fi
            else
                echo "No zip for \$gene / \$species; trying eutils fallback anyway..."
            fi
            download_and_normalize "\$species" "\$gene" "\$species_nospace"
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
                    echo "unzip failed for \$gene / \$species; trying eutils fallback anyway..."
                fi
            else
                echo "No zip for \$gene / \$species; trying eutils fallback anyway..."
            fi
            download_and_normalize "\$species" "\$gene" "\$species_nospace"
        done
    done

    # Concatenate all genes
    mito_genes_list=( ${params.mito_genes.collect { '"' + it + '"' }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { '"' + it + '"' }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )
    all_species_list=( ${params.species.collect { '"' + it + '"' }.join(' ')} )

    for gene in "\${all_genes_list[@]}"; do
        touch "data/concatenated_\${gene}.fasta"
        for species in "\${all_species_list[@]}"; do
            species_nospace=\$(echo "\$species" | tr ' ' '_')
            gene_file="data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna"
            if [ -s "\$gene_file" ]; then
                cat "\$gene_file" >> "data/concatenated_\${gene}.fasta"
            fi
        done
    done
    """
}

/*
 * Fallback downloader for queries missing from datasets download.
 * It searches NCBI nuccore and retrieves the first matching nucleotide entry with efetch.
 */
process nuccore_efetch_fallback {
    publishDir "${params.outDir}/fallback_${workflow.runName}", mode: 'copy', overwrite: true

    container 'thetinhornfish/efetch_retstart:24.8'

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

    # curl with retries + exponential backoff and jitter
    curl_retry() {
        max_attempts=3
        attempt=1
        while true; do
            curl -sG "\$@" && return 0
            if [ "\$attempt" -ge "\$max_attempts" ]; then
                return 1
            fi
            backoff=\$((2 ** (attempt - 1)))
            jitter=\$((RANDOM % 3))
            sleep_time=\$((backoff + jitter))
            sleep \$sleep_time
            attempt=\$((attempt + 1))
        done
    }

    efetch_try() {
        local out_file="\$1"
        url_base="\$2"
        seq_id="\$3"
        sleep \$((RANDOM % 2))

        # first try fasta_cds_na (preferred)
        if curl -fsSL -o "\${out_file}.tmp" "\${url_base}&id=\${seq_id}&rettype=fasta_cds_na&retmode=text"; then
            if [ -s "\${out_file}.tmp" ]; then
                mv "\${out_file}.tmp" "\${out_file}" && return 0
            else
                rm -f "\${out_file}.tmp"
            fi
        fi

        sleep \$((RANDOM % 2))
        # fallback to full fasta
        if curl -fsSL -o "\${out_file}.tmp" "\${url_base}&id=\${seq_id}&rettype=fasta&retmode=text"; then
            if [ -s "\${out_file}.tmp" ]; then
                mv "\${out_file}.tmp" "\${out_file}" && return 0
            else
                rm -f "\${out_file}.tmp"
            fi
        fi

        return 1
    }

    missing_file="data/missing_queries.tsv"
    touch data/efetch_recovered.tsv
    touch data/efetch_still_missing.tsv

    declare -A expected_len
    if [ -s data/genes_lengths_before.tsv ]; then
        while IFS=\$(printf '\\t') read -r _sp _gene _len; do
            [ -z "\$_gene" ] && continue
            [ -z "\$_len" ] && continue
            if [ -z "\${expected_len[\$_gene]:-}" ]; then
                expected_len[\$_gene]=\$_len
            fi
        done < data/genes_lengths_before.tsv
    fi

    if [ -s "\$missing_file" ]; then
        while IFS=\$(printf '\\t') read -r species_nospace gene; do
            [ -z "\$species_nospace" ] && continue
            [ -z "\$gene" ] && continue

            species="\${species_nospace//_/ }"
            out_dir="data/\$species_nospace/\${gene}/ncbi_dataset/data"
            out_file="\$out_dir/gene.fna"

            if [ -s "\$out_file" ]; then
                continue
            fi

            mkdir -p "\$out_dir"

            target_len="\${expected_len[\$gene]:-}"
            if [ -z "\$target_len" ]; then
                target_len=3000
            fi
            min_len=\$((\$target_len * 4 / 100))
            max_len=\$((\$target_len * 200 / 100))

            query="\${gene}[Gene Name] AND \${species}[Organism]"

            if ! search_xml=\$(curl_retry --data-urlencode "db=nuccore" --data-urlencode "retmax=20" --data-urlencode "term=\$query" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" 2>/dev/null); then
                search_xml=""
            fi

            id_list=\$(printf '%s' "\$search_xml" | grep -oE '<Id>[0-9]+</Id>' | sed 's:<Id>::g; s:</Id>::g')

            selected=0
            if [ -n "\$id_list" ]; then
                while IFS= read -r seq_id; do
                    [ -z "\$seq_id" ] && continue

                    tmp_fasta="\$out_file.tmp"
                    if ! efetch_try "\$tmp_fasta" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore" "\$seq_id"; then
                        rm -f "\$tmp_fasta" || true
                        continue
                    fi

                    # Filter to only the record matching [gene=GENE]
                    awk -v gene="\$gene" '
                        /^>/ { keep = (tolower(\$0) ~ tolower("\\[gene=" gene "\\]")) }
                        keep { print }
                    ' "\$tmp_fasta" > "\${tmp_fasta}.filtered"

                    if [ -s "\${tmp_fasta}.filtered" ]; then
                        mv "\${tmp_fasta}.filtered" "\$tmp_fasta"
                    else
                        rm -f "\$tmp_fasta" "\${tmp_fasta}.filtered"
                        continue
                    fi

                    seq_len=\$(awk '
                        BEGIN { n=0; seen=0 }
                        /^>/ { if (seen==1) exit; seen=1; next }
                        seen { gsub(/[[:space:]]/, "", \$0); n += length(\$0) }
                        END { print n }
                    ' "\$tmp_fasta")

                    if [ "\$seq_len" -lt "\$min_len" ] || [ "\$seq_len" -gt "\$max_len" ]; then
                        rm -f "\$tmp_fasta"
                        continue
                    fi

                    awk -v species="\$species" -v gene="\$gene" '
                        BEGIN { header_seen=0; done=0 }
                        /^>/ {
                            if (header_seen==1) { done=1; exit }
                            header_seen=1
                            print \$0 " [organism=" species "] [gene=" gene "]"
                            next
                        }
                        { if (!done) print }
                    ' "\$tmp_fasta" > "\$out_file"

                    rm -f "\$tmp_fasta"

                    if [ -s "\$out_file" ]; then
                        printf '%s\t%s\t%s\t%s\n' "\$species_nospace" "\$gene" "\$seq_id" "\$seq_len" >> data/efetch_recovered.tsv
                        selected=1
                        break
                    fi
                done <<< "\$id_list"
            fi

            # Third-level fallback: species-wide GenBank inspection
            if [ "\$selected" -eq 0 ]; then
                echo "No suitable gene-specific nuccore hit for \$species / \$gene; trying species-wide GenBank inspection..."

                if ! species_search_xml=\$(curl_retry --data-urlencode "db=nuccore" --data-urlencode "retmax=200" --data-urlencode "term=\${species}[Organism]" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" 2>/dev/null); then
                    species_search_xml=""
                fi

                species_id_list=\$(printf '%s' "\$species_search_xml" | grep -oE '<Id>[0-9]+</Id>' | sed 's:<Id>::g; s:</Id>::g')

                if [ -n "\$species_id_list" ]; then
                    while IFS= read -r seq_id; do
                        [ -z "\$seq_id" ] && continue

                        gb_tmp="\$out_file.\$seq_id.gb"
                        if ! curl_retry -o "\$gb_tmp" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\${seq_id}&rettype=gb&retmode=text" 2>/dev/null; then
                            rm -f "\$gb_tmp" || true
                            continue
                        fi

                        if grep -qiE "/gene=\"\${gene}\"|/product=.*\${gene}|\${gene}" "\$gb_tmp"; then
                            if ! efetch_try "\${out_file}.tmp" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore" "\${seq_id}"; then
                                rm -f "\${out_file}.tmp" "\$gb_tmp" || true
                                continue
                            fi

                            if [ -s "\${out_file}.tmp" ]; then
                                awk -v species="\$species" -v gene="\$gene" '
                                    BEGIN { header_seen=0; done=0 }
                                    /^>/ {
                                        if (header_seen==1) { done=1; exit }
                                        header_seen=1
                                        print \$0 " [organism=" species "] [gene=" gene "]"
                                        next
                                    }
                                    { if (!done) print }
                                ' "\${out_file}.tmp" > "\$out_file"

                                rm -f "\${out_file}.tmp" "\$gb_tmp"

                                if [ -s "\$out_file" ]; then
                                    seq_len=\$(awk '
                                        BEGIN { n=0; seen=0 }
                                        /^>/ { if (seen==1) exit; seen=1; next }
                                        seen { gsub(/[[:space:]]/, "", \$0); n += length(\$0) }
                                        END { print n }
                                    ' "\$out_file")

                                    if [ -n "\$seq_len" ] && [ "\$seq_len" -ge 100 ]; then
                                        printf '%s\t%s\t%s\t%s\n' "\$species_nospace" "\$gene" "\$seq_id" "\$seq_len" >> data/efetch_recovered.tsv
                                        selected=1
                                        break
                                    else
                                        rm -f "\$out_file"
                                    fi
                                fi
                            fi
                        fi
                        rm -f "\$gb_tmp"
                    done <<< "\$species_id_list"
                fi

                if [ "\$selected" -eq 0 ]; then
                    echo "No suitable nuccore hit for \$species / \$gene after species-wide search"
                    printf '%s\t%s\n' "\$species_nospace" "\$gene" >> data/efetch_still_missing.tsv
                    rm -f "\$out_file"
                fi
            fi
        done < "\$missing_file"
    fi

    # Rebuild concatenated FASTA files after fallback recovery.
    mito_genes_list=( ${params.mito_genes.collect { '"' + it + '"' }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { '"' + it + '"' }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )
    species_list=( ${params.species.collect { '"' + it + '"' }.join(' ')} )

    for gene in "\${all_genes_list[@]}"; do
        out_concat="data/concatenated_\${gene}.fasta"
        > "\$out_concat"
        for species in "\${species_list[@]}"; do
            species_nospace=\$(echo "\$species" | tr ' ' '_')
            gene_file="data/\$species_nospace/\${gene}/ncbi_dataset/data/gene.fna"
            if [ -s "\$gene_file" ]; then
                cat "\$gene_file" >> "\$out_concat"
            fi
        done
    done
    """
}

/*
 * This process removes duplicate sequences, keeping only the longest sequence for each species per gene
 */
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
    mito_genes_list=( ${params.mito_genes.collect { "'" + it + "'" }.join(' ')} )
    nuclear_genes_list=( ${params.nuclear_genes.collect { "'" + it + "'" }.join(' ')} )
    all_genes_list=( "\${mito_genes_list[@]}" "\${nuclear_genes_list[@]}" )
    species_list=( ${params.species.collect { '"' + it + '"' }.join(' ')} )
    expected_species_count=\${#species_list[@]}

    # Validation: Check that each gene has exactly one copy from all species
    echo "Validating gene sequence counts before alignment..."
    for gene in "\${all_genes_list[@]}"; do
        in_file="${inputDir}/concatenated_\${gene}.fasta"
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
    before_alignment(species=params.species)
    nuccore_efetch_fallback(before_alignment.out)
    deduplicate_sequences(nuccore_efetch_fallback.out)
    alignment(deduplicate_sequences.out)
    after_alignment(alignment.out)
    iqtree(after_alignment.out[0], after_alignment.out[3])
}