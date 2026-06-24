version 1.0

task array_to_zip {
    input {
        Array[File] input_files
        String zip_name = "compressed.zip"
    }

    command <<<
        echo "Compressing files into ~{zip_name}..."
        zip ~{zip_name} ~{sep=" " input_files}
    >>>

    output {
        File zipped_file = zip_name
    }

    runtime {
        docker: "python:3.11"
    }
}

task file_to_zip {
    input {
        File input_file
        String zip_name = "compressed.zip"
    }

    command <<<
        echo "Compressing ~{input_file} into ~{zip_name}..."
        zip ~{zip_name} ~{input_file}
    >>>

    output {
        File zipped_file = zip_name
    }

    runtime {
        docker: "python:3.11"
    }
}
