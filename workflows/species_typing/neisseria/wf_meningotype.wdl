version 1.0

import "../../../tasks/species_typing/neisseria/task_meningotype.wdl" as mt   # adapte le chemin

workflow wf_meningotype_fasta {
  meta {
    description: "N. meningitidis typing on genome assembly data with meningotype"
    author: "David Maimoun"
  }

  input {
    String samplename
    File assembly_fasta
  }

  call mt.meningotype as meningotype {
    input:
      assembly = assembly_fasta,
      samplename = samplename
  }

  output {
    File meningotype_tsv = meningotype.meningotype_tsv
    String meningotype_version = meningotype.meningotype_version

    String serogroup = meningotype.meningotype_serogroup
    String mlst = meningotype.meningotype_MLST
    String PorA = meningotype.meningotype_PorA
    String FetA = meningotype.meningotype_FetA
    String PorB = meningotype.meningotype_PorB
    String fHbp = meningotype.meningotype_fHbp
    String NHBA = meningotype.meningotype_NHBA
    String NadA = meningotype.meningotype_NadA
    String BAST = meningotype.meningotype_BAST
  }
}