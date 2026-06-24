version 1.0

workflow wf_unzip_into_bucket {
    input {
        File archive_file
        String output_bucket_path
        # example:
        # gs://my-bucket/results/run42
    }

    call unzip_and_rename {
        input:
            archive_file = archive_file
    }

    call upload_to_bucket {
        input:
            files = unzip_and_rename.renamed_files,
            output_bucket_path = output_bucket_path
    }
}

task unzip_and_rename {
    input {
        File archive_file
    }

    command <<<
        set -euo pipefail

        mkdir extracted
        mkdir renamed

        unzip "~{archive_file}" -d extracted

        find extracted -type f | while read f; do
            base=$(basename "$f")

            # remove everything before and including rerun_
            newname="${base#*rerun_}"

            cp "$f" "renamed/$newname"
        done
    >>>

    output {
        Array[File] renamed_files = glob("renamed/*")
    }

    runtime {
        docker: "ubuntu:22.04"
    }
}

task upload_to_bucket {
    input {
        Array[File] files
        String output_bucket_path
    }

    command <<<
        set -euo pipefail

        mkdir uploaded

        for f in ~{sep=' ' files}; do
            cp "$f" uploaded/
        done
    >>>

    output {
        Array[File] uploaded_files = glob("uploaded/*")
    }

    runtime {
        docker: "google/cloud-sdk:slim"
    }

    parameter_meta {
        uploaded_files: "Configure backend to delocalize to bucket"
    }
}