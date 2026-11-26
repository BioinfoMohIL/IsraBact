version 1.0

import "../../../tasks/species_typing/neisseria/task_neisseria_typing.wdl" as nm_typing

workflow wf_neisseria_typing_fasta {
  meta {
    description: "N. Meningitidis typing on genome assembly data"
    author: "David Maimoun"
  }

  input {
    String samplename
    File assembly_fasta 
    String pubmlst_api_key
  }
  
  call nm_typing.neisseria_typing as typing {
    input:
      assembly = assembly_fasta,
      pubmlst_api_key = pubmlst_api_key,
      samplename = samplename
  }

  call nm_typing.serogrouping as serogrouping {
      input: 
        assembly = assembly_fasta,
        samplename = samplename,
  }

  output {
      String neisseria_typing_version = typing.version
      String neisseria_typing_date = typing.date
      File nt_report = typing.typing_results
      File? nt_api_raw_results = typing.api_json_tar
      File? nt_logs = typing.logs


      String st = typing.st 
      String bast_type = typing.bast_type
      String clonal_complex = typing.clonal_complex
      String mendevar_bexsero_reactivity  = typing.bexsero_cross_reactivity      
      String mendevar_trumenba_reactivity = typing.trumenba_cross_reactivity
      String FetA_VR = typing.FetA_VR
      String fHbp_peptide = typing.fHbp_peptide
      String NHBA_peptide = typing.NHBA_peptide
      String NadA_peptide = typing.NadA_peptide
      String PorA_VR1 = typing.PorA_VR1 
      String PorA_VR2 = typing.PorA_VR2 
      String abcZ = typing.abcZ
      String adk = typing.adk
      String aroE = typing.aroE
      String gdh = typing.gdh
      String fumC = typing.fumC
      String pdhC = typing.pdhC
      String pgm = typing.pgm

      String serogroup = serogrouping.serogroup
      String genogroup = serogrouping.genogroup

      File? nt_loci_best_matches = typing.best_matches_tar

      String mlst_profile = typing.mlst_profile
      String bast_profile = typing.bast_profile
      String finetyping_profile = typing.finetyping_profile


  }
   
    
}


