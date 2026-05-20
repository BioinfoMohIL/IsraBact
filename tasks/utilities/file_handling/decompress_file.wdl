version 1.0

task decompress {
    input {
        File archive_file
        String output_dir = "decompressed"
    }

    command <<<
        set -e
        echo "Decompressing ~{archive_file} into ~{output_dir}..."
        mkdir -p "~{output_dir}"

        filename=$(basename "~{archive_file}")

        case "$filename" in
            *.zip)
                unzip -q "~{archive_file}" -d "~{output_dir}"
                ;;
            *.tar)
                tar -xf "~{archive_file}" -C "~{output_dir}"
                ;;
            *.tar.gz|*.tgz)
                tar -xzf "~{archive_file}" -C "~{output_dir}"
                ;;
            *.tar.bz2)
                tar -xjf "~{archive_file}" -C "~{output_dir}"
                ;;
            *.gz)
                # Cas particulier : extrait un SEUL fichier directement à la racine du dossier
                # On garde le sous-dossier pour éviter les conflits de scope Cromwell
                gzip -d -c "~{archive_file}" > "~{output_dir}/$(basename "$filename" .gz)"
                ;;
            *)
                echo "Unsupported file type: $filename"
                exit 1
                ;;
        esac
    >>>

    output {
        String decompressed_dir_path = output_dir
        Array[File] decompressed_files = glob("~{output_dir}/*")
        
        # On prend directement le premier élément de l'array non-optionnel
        File decompressed_single_file = glob("~{output_dir}/*")[0]
    }

    runtime {
        docker: "python:3.11"
        cpu: 1
        memory: "2 GB"
    }
}