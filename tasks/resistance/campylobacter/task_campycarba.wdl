version 1.0

task campycarba {
  input {
    File assembly_fasta
    String samplename
    String species                    # must contain "jejuni" or "coli" (selects the default reference set in the image)

    File? porA_reference_faa          # optional override: MATURE MOMP protein (position 1 = first residue after the signal peptide)
    File? oxa_reference_fasta         # optional override: blaOXA-61 locus (upstream + CDS); header must contain cds_start=N
    File? knowledge_base_tsv          # defaults to the copy baked into the image
    File? bam                         # optional: coordinate-sorted BAM mapped to assembly_fasta (read-depth duplication check)
    String docker = "bioinfomoh/campycarba:1.1"

    Int cpu = 4
    Int memory_gb = 4
    Int disk_size_gb = 20
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    campycarba --help > /dev/null
    mkdir -p out
    campycarba \
      --assembly "~{assembly_fasta}" \
      --sample "~{samplename}" \
      --species "~{species}" \
      --kb "~{if defined(knowledge_base_tsv) then select_first([knowledge_base_tsv]) else '/opt/campycarba/campycarba_kb.tsv'}" \
      ~{if defined(porA_reference_faa) then "--porA-ref '" + select_first([porA_reference_faa]) + "'" else ""} \
      ~{if defined(oxa_reference_fasta) then "--oxa-ref '" + select_first([oxa_reference_fasta]) + "'" else ""} \
      ~{if defined(bam) then "--bam '" + select_first([bam]) + "'" else ""} \
      --outdir out
  >>>

  output {
    File summary_tsv  = "out/~{samplename}_summary.tsv"
    File variants_tsv = "out/~{samplename}_variants.tsv"
    File report_txt   = "out/~{samplename}_report.txt"
    String status     = read_string("out/~{samplename}_status.txt")
    # Alignments / translated ORF for a visual check in MEGA, AliView or Jalview (empty if the gene was not found)
    Array[File] alignment_files = flatten([glob("out/*_alignment.fasta"), glob("out/*_porA_protein.faa")])
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} HDD"
    preemptible: preemptible
  }
}
