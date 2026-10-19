version 1.0

import "../../../tasks/resistance/campylobacter/task_campycarba.wdl" as task_campycarba

workflow wf_campycarba {
  input {
    File assembly_fasta
    String samplename
    String species                    # e.g. "Campylobacter jejuni" or "Campylobacter coli"

    File? porA_reference_faa
    File? oxa_reference_fasta
    File? knowledge_base_tsv
    File? bam
    String docker = "bioinfomoh/campycarba:1.1"
  }

  call task_campycarba.campycarba {
    input:
      assembly_fasta = assembly_fasta,
      samplename = samplename,
      species = species,
      porA_reference_faa = porA_reference_faa,
      oxa_reference_fasta = oxa_reference_fasta,
      knowledge_base_tsv = knowledge_base_tsv,
      bam = bam,
      docker = docker
  }

  output {
    File campycarba_summary_tsv  = campycarba.summary_tsv
    File campycarba_variants_tsv = campycarba.variants_tsv
    File campycarba_report_txt   = campycarba.report_txt
    String campycarba_status     = campycarba.status
    Array[File] campycarba_alignment_files = campycarba.alignment_files
  }
}
