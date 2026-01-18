version 1.0

task decompress_file {
    input {
        File archive_file
        String output_dir = "decompressed"
    }

    command <<<
        echo "Decompressing ~{archive_file} into ~{output_dir}..."
        mkdir -p ~{output_dir}

        filename=$(basename ~{archive_file})

        case "$filename" in
            *.zip)
                unzip ~{archive_file} -d ~{output_dir}
                ;;
            *.tar)
                tar -xf ~{archive_file} -C ~{output_dir}
                ;;
            *.tar.gz|*.tgz)
                tar -xzf ~{archive_file} -C ~{output_dir}
                ;;
            *.tar.bz2)
                tar -xjf ~{archive_file} -C ~{output_dir}
                ;;
            *.gz)
                gzip -d -c ~{archive_file} > ~{output_dir}/$(basename "$filename" .gz)
                ;;
            *)
                echo "Unsupported file type: $filename"
                exit 1
                ;;
        esac
    >>>


    output {
        String decompressed_dir = output_dir
    }

    runtime {
        docker: "python:3.11"
    }
}
