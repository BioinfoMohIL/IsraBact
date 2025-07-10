version 1.0


task allele_calling {
    input {
        File schema_zip
        File prodigal_zip
        Array[File]? assemblies
        File? assemblies_zipped
        Int cpu = 20
    }

    command <<<
        echo "Decompressing files..."

        for archive_file in ~{schema_zip} ~{prodigal_zip}; do
            output_dir=$(basename "$archive_file")_unzipped
            echo "Decompressing $archive_file into $output_dir"
            mkdir -p "$output_dir"

            filename=$(basename "$archive_file")

            case "$filename" in
                *.zip)
                    unzip "$archive_file" -d "$output_dir"
                    ;;
                *.tar)
                    tar -xf "$archive_file" -C "$output_dir"
                    ;;
                *.tar.gz|*.tgz)
                    tar -xzf "$archive_file" -C "$output_dir"
                    ;;
                *.tar.bz2)
                    tar -xjf "$archive_file" -C "$output_dir"
                    ;;
                *.gz)
                    gzip -d -c "$archive_file" > "$output_dir"/$(basename "$filename" .gz)
                    ;;
                *)
                    echo "Unsupported file type: $filename"
                    exit 1
                    ;;
            esac
        done

        echo "[1] Adapt your external schema"
        chewBBACA.py PrepExternalSchema \
                    -i $(basename ~{schema_zip})_unzipped \
                    -o schema_adapted \
                    --ptf $(basename ~{prodigal_zip})_unzipped/*.trn \
                    --cpu ~{cpu}  

        echo "[2] Preparing assemblies input"

        mkdir -p assemblies_files

        # Case 1 : assemblies_zipped 
        if [ -s "~{assemblies_zipped}" ]; then
            echo "Decompressing assemblies_zipped"

            tmp_dir=$(mktemp -d)

            # Try unzip
            unzip "~{assemblies_zipped}" -d "$tmp_dir" || \
            tar -xf "~{assemblies_zipped}" -C "$tmp_dir" || \
            tar -xzf "~{assemblies_zipped}" -C "$tmp_dir" || \
            { echo "Unknown compression format for assemblies_zipped"; exit 1; }

            echo "Moving files to assemblies_files"
            find "$tmp_dir" -type f -exec mv {} assemblies_files/ \;

            rm -rf "$tmp_dir"
        fi
        
        # Case 2 : assemblies array files
        if [ "~{sep=' ' assemblies}" != "" ]; then
            echo "Copying assemblies list"
            for f in ~{sep=' ' assemblies}; do
                cp "$f" assemblies_files/
            done
        fi

        echo "[3] Running AlleleCall"
        chewBBACA.py AlleleCall -i assemblies_files -g schema_adapted -o calling --cpu ~{cpu}

        echo "[4] Cleaning sample names"
        results_file=$(find calling -name "results_alleles.tsv" | head -n 1)
        
        awk -F'\t' 'BEGIN{OFS="\t"} {
            name = $1
            sub(/\.[^.]+$/, "", name)
            split(name, parts, /[_-]/)
            $1 = parts[1]
            print
        }' "$results_file" > alleles_cleaned.tsv

        echo "Compressing output ..."
        tar -czf schema_adapted.tar.gz schema_adapted
        tar -czf calling.tar.gz calling

    >>>

    output {
        File calling_dir = "calling.tar.gz"
        File alleles_cleaned = "alleles_cleaned.tsv"
        File schema_adapted = "schema_adapted.tar.gz"
    }

    runtime {
        docker: "bioinfomoh/run_chewbbaca:1"
        cpu : "~{cpu}"
    }
}

task extract_cgmlst {
    input {
        File cleaned_results
    }

    command {
        echo "[4] Extraction for data visualization"
        chewBBACA.py ExtractCgMLST -i ~{cleaned_results} -o visualization

        tar -czf visualization.tar.gz visualization

    }

    output {
        File visualization_zip = "visualization.tar.gz"
        File visualization_file = "visualization/cgMLST.tsv"
    }

    runtime {
        docker: "bioinfomoh/run_chewbbaca:1"
    }
}
