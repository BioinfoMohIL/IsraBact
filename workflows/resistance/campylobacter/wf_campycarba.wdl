version 1.0

import "../../../tasks/resistance/campylobacter/task_campycarba.wdl" as task_campycarba
import "../../../tasks/resistance/campylobacter/task_campycarba_pair.wdl" as task_campycarba_pair

workflow wf_campycarba {
  input {
    File assembly_fasta
    String samplename
    String species                    # e.g. "Campylobacter jejuni" or "Campylobacter coli"

    # Optional: give the SUSCEPTIBLE isolate of the same patient/lineage to also run the paired comparison
    File? susceptible_assembly_fasta
    String? susceptible_samplename

    File? porA_reference_faa
    File? oxa_reference_fasta
    File? knowledge_base_tsv
    File? bam
    String docker = "bioinfomoh/campycarba:1.0"
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

  if (defined(susceptible_assembly_fasta)) {
    call task_campycarba_pair.campycarba_pair {
      input:
        resistant_assembly_fasta = assembly_fasta,
        samplename = samplename,
        susceptible_assembly_fasta = select_first([susceptible_assembly_fasta]),
        susceptible_samplename = select_first([susceptible_samplename, "susceptible"]),
        species = species,
        porA_reference_faa = porA_reference_faa,
        oxa_reference_fasta = oxa_reference_fasta,
        docker = docker
    }
  }

  output {
    File campycarba_summary_tsv  = campycarba.summary_tsv
    File campycarba_variants_tsv = campycarba.variants_tsv
    File campycarba_report_txt   = campycarba.report_txt
    String campycarba_status     = campycarba.status
    Array[File] campycarba_alignment_files = campycarba.alignment_files

    # Only present when a susceptible isolate was given
    File? campycarba_pair_report_txt   = campycarba_pair.pair_report_txt
    File? campycarba_pair_variants_tsv = campycarba_pair.pair_variants_tsv
    Array[File]? campycarba_pair_alignment_files = campycarba_pair.pair_alignment_files
  }
}
