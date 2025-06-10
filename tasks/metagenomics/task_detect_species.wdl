version 1.0


task detect_species {
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
            kraken2 $mode $compressed --threads "~{cpu}" --use-names --db /app/db/kraken_db \
                --report "~{sample_id}.report" --paired "~{read1}" "~{read2}" --output -

            echo "~{sample_id}" > ~{sample_id}_sample_detected.csv
            
            declare -A species
            species["NM"]="Neisseria Meningitidis"
            species["NG"]="Neisseria Gonorrhoeae"
            species["HI"]="Haemophilus Influenzae"
            species["SH"]="Salmonella"
            species["SO"]="Salmonella"
            species["LC"]="Listeria monocytogenes"
            species["LF"]="Listeria monocytogenes"
            species["SG"]="Shigella"
            species["CA"]="Campylobacter"
            species["VIB"]="Vibrio"
            species["V"]="Vibrio"
            species["EC"]="Escherichia coli"
            species["SA"]="Staphylococcus aureus"
            species["BP"]="Bordetella pertussis"
            species["SP"]="Streptococcus pneumoniae"
            species["ST"]="Streptococcus pyogenes"
            species["ST"]="Streptococcus agalactiae"
            species["LG"]="Legionella pneumophila"
            species["LW"]="Legionella pneumophila"
            species["CB"]="Corynebacterium diphtheriae"
            species["HI"]="Haemophilus influenzae"
            species["NM"]="Neisseria meningitidis"
            species["M" ]="Neisseria meningitidis" 

            prefix=$(echo "~{sample_id}" | grep -o '^[^0-9]*')

            # Extract detected species from report
            detected=$(awk -F'\t' '$4 == "S" {gsub(/^[ \t]+/, "", $6); print $6; exit}' "~{sample_id}.report")

            # Convert both to lowercase for case-insensitive comparison
            detected_lower=$(echo "$detected" | tr '[:upper:]' '[:lower:]')
            expected_lower=$(echo "${species[$prefix]}" | tr '[:upper:]' '[:lower:]')

            if [[ "$prefix" == "ST" && "$detected" == *"Streptococcus"* ]]; then
                # handle ST samples (can be pyogenes or agalactiae)
                echo "Met ST , and detected Streptococcus"
                echo "~{sample_id},${detected},+" > ~{sample_id}_sample_detected.csv

            elif [[ "$prefix" == "SG" && "$detected" == "Escherichia coli" ]]; then
                # handle Shigella (SG) - Kraken detect it as EC
                echo "~{sample_id},${detected},??" > ~{sample_id}_sample_detected.csv

            elif [[ "$detected_lower" == *"$expected_lower"* ]]; then
                echo "~{sample_id},${detected},+" > ~{sample_id}_sample_detected.csv

            else   
                if [[ -n "${detected}" && "${detected}" != "" ]]; then
                    echo "~{sample_id},${detected},xxx" > ~{sample_id}_sample_detected.csv
                else
                    echo "~{sample_id},,xxx" > ~{sample_id}_sample_detected.csv
                fi
            fi
        fi
    >>>

    output {
        File report = "~{sample_id}.report"
        String sample_detected = read_string("~{sample_id}_sample_detected.csv")
    }

    runtime {
        docker: docker
        cpu: cpu
    }
}