version 1.0


task kraken2 {
  input {
    File read1
    File  read2
    String sample_id
    File kraken_db = "gs://theiagen-public-resources-rp/reference_data/databases/kraken2/k2_viral-refseq_human-GRCh38_20260220.tar.gz"
    
    # String docker = "bioinfomoh/specie_detection:1"
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/kraken2:2.17.1"

    Int cpu = 20
  
  }

  command <<<
        mode=""
        compressed=""

        # Check if paired mode should be used
        if ! [ -z "~{read2}" ]; then
            echo "Reads are paired..."
            mode="--paired"
        fi

        # Determine if reads are compressed
        if [[ "~{read1}" == *.gz ]]; then
            echo "Reads are compressed..."
            compressed="--gzip-compressed"
        fi

        tar -xzf "~{kraken_db}" -C kraken_db

        echo "DB Contents:"
        find kraken_db -maxdepth 2


        # Run Kraken2
        echo "Running Kraken2..."

        kraken2 -v | awk '/Kraken/ {print "Kraken v" $3}' | tee VERSION
        # /app/db/kraken_db
        kraken2 $mode $compressed --threads "~{cpu}" --use-names --db kraken_db \
            --report "~{sample_id}_report.txt" --paired "~{read1}" "~{read2}" --output -

        detected=$(awk -F'\t' '$4 == "S" {gsub(/^[ \t]+/, "", $6); print $6; exit}' "~{sample_id}_report.txt")

        echo "${detected}" > ~{sample_id}_detected.txt    
    >>>

    output {
        String version = read_string("VERSION")
        File report = "~{sample_id}_report.txt"
        String species_detected = read_string("~{sample_id}_detected.txt")
    }

    runtime {
        docker: docker
        cpu: cpu
    }
}

task detect_species_bs_reads {
  input {
    File? read1
    File? read2
    String sample_id
    
    String docker = "bioinfomoh/specie_detection:1"
    Int cpu = 20
    
  }

  command <<<
        mode=""
        compressed=""

        if [ ! -s "~{read1}" ]; then
            echo "~{sample_id} not found" > "~{sample_id}.report"
            echo "~{sample_id},not_found,xxx" > species_detected.csv        
        else
            # Check if paired mode should be used
            if ! [ -z "~{read2}" ]; then
                echo "Reads are paired..."
                mode="--paired"
            fi

            # Determine if reads are compressed
            if [[ "~{read1}" == *.gz ]]; then
                echo "Reads are compressed..."
                compressed="--gzip-compressed"
            fi

            # Run Kraken2
            echo "Running Kraken2..."

            kraken2 -v | grep -oP '^Kraken\s+\Kversion\s+\K[0-9.]+' | awk '{print "Kraken v"$1}' | tee VERSION

            kraken2 $mode $compressed --threads "~{cpu}" --use-names --db /app/db/kraken_db \
                --report "~{sample_id}_report.txt" --paired "~{read1}" "~{read2}" --output -

            detected=$(awk -F'\t' '$4 == "S" {gsub(/^[ \t]+/, "", $6); print $6; exit}' "~{sample_id}_report.txt")

            echo "Sample,Detected,Match" > ~{sample_id}_detected.csv
            echo "~{sample_id}","${detected}" >> ~{sample_id}_detected.csv
        fi
          
    >>>

    output {
        String version = read_string("VERSION")
        File? report = "~{sample_id}_report.txt"
        File sample_detected = "~{sample_id}_detected.csv"
    }

    runtime {
        docker: docker
        cpu: cpu
    }
}