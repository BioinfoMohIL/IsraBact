version 1.0

task campycarba_pair {
  input {
    File resistant_assembly_fasta
    String samplename                          # the RESISTANT isolate
    File susceptible_assembly_fasta
    String susceptible_samplename
    String species                             # must contain "jejuni" or "coli" (selects the default reference set)

    File? porA_reference_faa                   # optional overrides (same files as for campycarba)
    File? oxa_reference_fasta
    String docker = "bioinfomoh/campycarba:1.2"

    Int cpu = 2
    Int memory_gb = 4
    Int disk_size_gb = 20
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    mkdir -p out
    compare_pair \
      --susceptible "~{susceptible_assembly_fasta}" \
      --resistant "~{resistant_assembly_fasta}" \
      --name-s "~{susceptible_samplename}" \
      --name-r "~{samplename}" \
      --species "~{species}" \
      ~{if defined(porA_reference_faa) then "--porA-ref '" + select_first([porA_reference_faa]) + "'" else ""} \
      ~{if defined(oxa_reference_fasta) then "--oxa-ref '" + select_first([oxa_reference_fasta]) + "'" else ""} \
      --prefix "~{samplename}_vs_~{susceptible_samplename}" \
      --outdir out
  >>>

  output {
    File pair_report_txt   = "out/~{samplename}_vs_~{susceptible_samplename}_pair_report.txt"
    File pair_variants_tsv = "out/~{samplename}_vs_~{susceptible_samplename}_pair_variants.tsv"
    Array[File] pair_alignment_files = glob("out/*_alignment.fasta")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} HDD"
    preemptible: preemptible
  }
}
