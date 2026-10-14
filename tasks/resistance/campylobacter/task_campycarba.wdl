version 1.0

task campycarba {
  input {
    File assembly_fasta
    String sample_name

    File? porA_reference_faa          # MATURE MOMP protein (position 1 = first residue after the signal peptide)
    File? oxa_reference_fasta         # blaOXA-61 locus: upstream region + CDS; header must contain cds_start=N
    File? knowledge_base_tsv         # defaults to the copy baked into the image
    String species = "unknown"
    File? bam                        # optional: coordinate-sorted BAM mapped to assembly_fasta (read-depth duplication check)
    String docker = "bioinfomoh/campycarba:1.0"
    
    Int cpu = 4
    Int memory_gb = 4
    Int disk_size_gb = 20
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    mkdir -p out
    campycarba \
      --assembly "~{assembly_fasta}" \
      --sample   "~{sample_name}" \
      --porA-ref "~{if defined(porA_reference_faa) then select_first([porA_reference_faa]) else '/opt/campycarba/porA_mature_ref.faa'}" \
      --oxa-ref  "~{if defined(oxa_reference_fasta) then select_first([oxa_reference_fasta]) else '/opt/campycarba/blaOXA-61_locus_ref.fasta'}" \
      --kb       "~{if defined(knowledge_base_tsv) then select_first([knowledge_base_tsv]) else '/opt/campycarba/campycarba_kb.tsv'}" \
      --species  "~{species}" \
      ~{"--bam " + bam} \
      --outdir out
  >>>

  output {
    File summary_tsv  = "out/~{sample_name}_summary.tsv"
    File variants_tsv = "out/~{sample_name}_variants.tsv"
    File report_txt   = "out/~{sample_name}_report.txt"
    String status     = read_string("out/~{sample_name}_status.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} HDD"
    preemptible: preemptible
  }
}
