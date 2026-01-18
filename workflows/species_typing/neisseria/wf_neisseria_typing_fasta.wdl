version 1.0

# import "../../../tasks/species_typing/neisseria/neisseria_typing.wdl" as nm_typing
import "../../../tasks/species_typing/neisseria/nt.wdl" as nm_typing

workflow wf_neisseria_typing_fasta {
  meta {
    description: "N. Meningitidis typing on genome assembly data"
    author: "David Maimoun"
  }

  input {
    String samplename
    File assembly_fasta 
  }
  
  call nm_typing.neisseria_typing as typing {
    input:
      assembly = assembly_fasta,
      samplename = samplename
  }

  call nm_typing.serogrouping as serogrouping {
      input: 
        assembly = assembly_fasta,
        samplename = samplename,
  }

  output {
      File typing_results = typing.typing_results
      String neisseria_typing_version = typing.version
      String st = typing.st 
      String bast_type = typing.bast_type
      String clonal_complex = typing.clonal_complex
      String bexsero_cross_reactivity  = typing.bexsero_cross_reactivity      
      String trumenba_cross_reactivity = typing.trumenba_cross_reactivity
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


  }
   
    
}


