process download_fallback {
    publishDir "${params.outDir}", mode: 'copy'

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

        if curl -fsSL -o "\${out_file}.tmp" "\${url_base}&id=\${seq_id}&rettype=fasta_cds_na&retmode=text"; then
            if [ -s "\${out_file}.tmp" ]; then
                mv "\${out_file}.tmp" "\${out_file}" && return 0
            else
                rm -f "\${out_file}.tmp"
            fi
        fi

        sleep \$((RANDOM % 2))
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
        while IFS=\$(printf '\t') read -r _sp _gene _len; do
            [ -z "\$_gene" ] && continue
            [ -z "\$_len" ] && continue
            if [ -z "\${expected_len[\$_gene]:-}" ]; then
                expected_len[\$_gene]=\$_len
            fi
        done < data/genes_lengths_before.tsv
    fi

    if [ -s "\$missing_file" ]; then
        while IFS=\$(printf '\t') read -r species_nospace gene; do
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
    """
}
