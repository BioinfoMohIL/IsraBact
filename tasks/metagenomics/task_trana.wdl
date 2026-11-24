version 1.0

task task_trana {

  input {
    File sample_sheet
    String outdir
    String db_path
    String seqtype
    String fastq_pass
    String? barcodes_sheet
  }

  command <<<

    set -euo pipefail

    echo "Starting TRANA pipeline"

    if [ -z "~{barcodes_sheet}" ]
    then
      nextflow run main.nf \
        --input ~{sample_sheet} \
        --outdir ~{outdir} \
        --db ~{db_path} \
        --seqtype ~{seqtype} \
        --merge_fastq_pass ~{fastq_pass} \
        -profile docker
    else
      nextflow run main.nf \
        --input ~{sample_sheet} \
        --outdir ~{outdir} \
        --db ~{db_path} \
        --seqtype ~{seqtype} \
        --merge_fastq_pass ~{fastq_pass} \
        --barcodes_samplesheet ~{barcodes_sheet} \
        -profile docker
    fi

  >>>

  output {
    Directory results = "~{outdir}"
  }

  runtime {
    docker: "quay.io/biocontainers/nextflow:latest"
  }

}