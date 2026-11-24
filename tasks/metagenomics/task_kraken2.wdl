version 1.0

task kraken2 {
  input {
    File read1
    File read2
    String sample_id
    File kraken_db = "gs://fc-5d4556f8-3de6-4709-85da-11445772644d/db/minikraken2_v2_8GB_201904_UPDATE.tar.gz"

    Int memory = 32
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/kraken2:2.17.1"
    Int cpu = 20
  }

  command <<<
        set -euo pipefail

        mode=""
        compressed=""

        if ! [ -z "~{read2}" ]; then
            echo "Reads are paired..."
            mode="--paired"
        fi

        if [[ "~{read1}" == *.gz ]]; then
            echo "Reads are compressed..."
            compressed="--gzip-compressed"
        fi

        if [ ! -s "~{kraken_db}" ]; then
            echo "ERROR: DB file missing or empty"
            exit 1
        fi

        mkdir -p kraken_db
        tar -xzf "~{kraken_db}" -C kraken_db
        find kraken_db -maxdepth 2

        echo "Running Kraken2..."
        kraken2 -v | awk '/Kraken/ {print "Kraken v" $3}' | tee VERSION

        kraken2 $mode $compressed --threads "~{cpu}" --use-names \
            --db kraken_db/minikraken2_v2_8GB_201904_UPDATE \
            --report "~{sample_id}_report.txt" "~{read1}" "~{read2}" --output -

        REPORT="~{sample_id}_report.txt"

        # --- Top genus / top species (highest %), format: name<TAB>pct ---
        top_taxon() {
            awk -F'\t' -v rank="$1" '
                $4 == rank {
                    name = $6; gsub(/^[ \t]+/, "", name)
                    pct = $1; gsub(/^[ \t]+/, "", pct)
                    if (pct + 0 > best + 0 || found == 0) { best = pct; bname = name; found = 1 }
                }
                END { if (found) printf "%s\t%s", bname, best }
            ' "$REPORT"
        }

        IFS=$'\t' read -r genus genus_pct <<< "$(top_taxon G)" || true
        IFS=$'\t' read -r species species_pct <<< "$(top_taxon S)" || true

        # Ex: "Bordetella:98.73, Bordetella pertussis:34.99"
        echo "${genus:-NA}:${genus_pct:-0}, ${species:-NA}:${species_pct:-0}" > top_taxa.txt

        # Kept for backward compatibility
        echo "${species}" > "~{sample_id}_detected.txt"

        # Genus alone, used by task_species_match (kraken2 and ML share the same match logic)
        echo "${genus}" > genus.txt


    >>>

    output {
        String version = read_string("VERSION")
        File report = "~{sample_id}_report.txt"
        String species_detected = read_string("~{sample_id}_detected.txt")
        String top_taxa = read_string("top_taxa.txt")
        String genus_detected = read_string("genus.txt")
    }

    runtime {
        docker: "~{docker}"
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB" # TES
        preemptible: 0
    }
}