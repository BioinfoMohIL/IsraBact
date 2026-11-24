version 1.0

import "../../../tasks/species_typing/neisseria/neisseria_typing.wdl" as nm_typing
import "../../../tasks/species_typing/neisseria/task_meningotype.wdl" as mt   # adapte le chemin

workflow wf_neisseria_typing_fasta {
  meta {
    description: "N. Meningitidis typing on genome assembly data (meningotype by default, PubMLST API optional)"
    author: "David Maimoun"
  }

  input {
    String samplename
    File assembly_fasta
    Boolean call_api = false
    String pubmlst_api_key = ""
  }

  Boolean use_api = call_api && (pubmlst_api_key != "")

  # ---------- PubMLST API ----------
  if (use_api) {
    call nm_typing.neisseria_typing as typing {
      input:
        assembly = assembly_fasta,
    }

    call nm_typing.serogrouping as serogrouping {
      input:
        assembly = assembly_fasta,
        samplename = samplename
    }
  }

  # ---------- Meningotype (default) ----------
  if (!use_api) {
    call mt.meningotype as mening {
      input:
        assembly = assembly_fasta,
        samplename = samplename
    }
  }

  output {
    # ===== API (PubMLST) =====
    String? nt_api_version = typing.version
    String? nt_api_date = typing.date

    String? nt_api_fHbp_peptide = typing.fHbp_peptide
    String? nt_api_NHBA_peptide = typing.NHBA_peptide
    String? nt_api_NadA_peptide = typing.NadA_peptide
    String? nt_api_PorA_VR1 = typing.PorA_VR1
    String? nt_api_PorA_VR2 = typing.PorA_VR2
    String? nt_api_abcZ = typing.abcZ
    String? nt_api_adk = typing.adk
    String? nt_api_aroE = typing.aroE
    String? nt_api_gdh = typing.gdh
    String? nt_api_fumC = typing.fumC
    String? nt_api_pdhC = typing.pdhC
    String? nt_api_pgm = typing.pgm
    String? nt_api_st = typing.st
    String? nt_api_FetA_VR = typing.FetA_VR
    String? nt_api_bast_type = typing.bast_type
    String? nt_api_clonal_complex = typing.clonal_complex
    String? nt_api_mendevar_bexsero_reactivity = typing.bexsero_cross_reactivity
    String? nt_api_mendevar_trumenba_reactivity = typing.trumenba_cross_reactivity
    String? nt_api_mlst_profile = typing.mlst_profile
    String? nt_api_bast_profile = typing.bast_profile
    String? nt_api_finetyping_profile = typing.finetyping_profile
    File? nt_api_loci_best_matches = typing.best_matches_tar

    String? nt_api_serogroup = serogrouping.serogroup
    String? nt_api_genogroup = serogrouping.genogroup

    # ===== Meningotype =====
    File? nt_mening_tsv = mening.meningotype_tsv
    String? nt_mening_version = mening.meningotype_version
    String? nt_mening_serogroup = mening.meningotype_serogroup
    String? nt_mening_mlst = mening.meningotype_MLST
    String? nt_mening_PorA = mening.meningotype_PorA
    String? nt_mening_FetA = mening.meningotype_FetA
    String? nt_mening_PorB = mening.meningotype_PorB
    String? nt_mening_fHbp = mening.meningotype_fHbp
    String? nt_mening_NHBA = mening.meningotype_NHBA
    String? nt_mening_NadA = mening.meningotype_NadA
    String? nt_mening_BAST = mening.meningotype_BAST
    String? nt_mening_bexsero_mendevar = mening.meningotype_Bexsero_MenDeVAR
    String? nt_mening_trumenba_mendevar = mening.meningotype_Trumenba_MenDeVAR
  }
}