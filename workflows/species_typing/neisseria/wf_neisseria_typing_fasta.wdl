version 1.0

import "../../../tasks/species_typing/neisseria/neisseria_typing.wdl" as nm_typing
import "../../../tasks/species_typing/neisseria/task_meningotype.wdl" as mt   # adapte le chemin

workflow wf_neisseria_typing_fasta {
  meta {
    description: "N. meningitidis typing on assembly: meningotype by default, PubMLST API optional"
    author: "David Maimoun"
  }

  input {
    String samplename
    File assembly_fasta
    Boolean call_api = false
    String pubmlst_api_key = ""
  }

  Boolean use_api = call_api && (pubmlst_api_key != "")

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

  if (!use_api) {
    call mt.meningotype as mening {
      input:
        assembly = assembly_fasta,
        samplename = samplename
    }
  }

  output {
    String nt_source = if use_api then "pubmlst_api" else "meningotype"

    # ---- Champs communs aux deux branches ----
    String nt_version   = select_first([typing.version, mening.meningotype_version])
    String nt_serogroup = select_first([serogrouping.serogroup, mening.meningotype_serogroup])
    String nt_st        = select_first([typing.st, mening.meningotype_MLST])
    String nt_fHbp      = select_first([typing.fHbp_peptide, mening.meningotype_fHbp])
    String nt_NHBA      = select_first([typing.NHBA_peptide, mening.meningotype_NHBA])
    String nt_NadA      = select_first([typing.NadA_peptide, mening.meningotype_NadA])
    String nt_FetA      = select_first([typing.FetA_VR, mening.meningotype_FetA])
    String nt_bast      = select_first([typing.bast_type, mening.meningotype_BAST])
    String nt_PorA      = if use_api
                          then select_first([typing.PorA_VR1, ""]) + "," + select_first([typing.PorA_VR2, ""])
                          else select_first([mening.meningotype_PorA, ""])

    # ---- Champs propres à une branche ----
    File?   nt_tsv              = mening.meningotype_tsv
    String? nt_PorB             = mening.meningotype_PorB

    String? nt_date             = typing.date
    String? nt_genogroup        = serogrouping.genogroup
    String? nt_clonal_complex   = typing.clonal_complex
    String? nt_abcZ             = typing.abcZ
    String? nt_adk              = typing.adk
    String? nt_aroE             = typing.aroE
    String? nt_gdh              = typing.gdh
    String? nt_fumC             = typing.fumC
    String? nt_pdhC             = typing.pdhC
    String? nt_pgm              = typing.pgm
    String? nt_bexsero_mendevar  = typing.bexsero_cross_reactivity
    String? nt_trumenba_mendevar = typing.trumenba_cross_reactivity
    String? nt_mlst_profile       = typing.mlst_profile
    String? nt_bast_profile       = typing.bast_profile
    String? nt_finetyping_profile = typing.finetyping_profile
    File?   nt_loci_best_matches  = typing.best_matches_tar
  }
}