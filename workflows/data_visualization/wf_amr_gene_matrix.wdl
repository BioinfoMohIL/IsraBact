version 1.0

import "../../../tasks/data_visualization/task_build_gene_matrix.wdl" as tasks

## wf_amr_gene_matrix.wdl
##
## 4 inputs optionnels, tous independants les uns des autres :
##   - abricate_reports             : rapports AbriCate (genes de virulence)
##   - amrfinder_amr_reports        : rapports AMRFinderPlus - genes de resistance (amr_report)
##   - amrfinder_stress_reports     : rapports AMRFinderPlus - genes de stress (stress_report)
##   - amrfinder_virulence_reports  : rapports AMRFinderPlus - genes de virulence (virulence_report)
##
## Seuls les inputs fournis sont traites : si un input est omis, la task correspondante
## n'est pas executee et ses outputs restent vides. Cela permet par exemple de ne
## generer qu'une matrice de virulence sans avoir a fournir les 3 autres.
##
## Pour chaque input fourni, produit :
##   - matrice YES/NO   (TSV + XLSX)
##   - matrice binaire  (TSV + XLSX)
##   - un fichier concatenant les rapports bruts d'entree (avec colonne Key)

workflow wf_amr_gene_matrix {
  input {
    Array[File]? abricate_hits_reports
    Array[File]? amrfinder_amr_reports
    Array[File]? amrfinder_stress_reports
    Array[File]? amrfinder_virulence_reports
  }

  if (defined(abricate_reports)) {
    call tasks.build_gene_matrix as abricate_matrix {
      input:
        tsv_reports = select_first([abricate_reports]),
        gene_column = "GENE",
        source_name = "abricate"
    }
  }

  if (defined(amrfinder_amr_reports)) {
    call tasks.build_gene_matrix as amrfinder_amr_matrix {
      input:
        tsv_reports = select_first([amrfinder_amr_reports]),
        gene_column = "Element symbol",
        source_name = "amrfinder_amr"
    }
  }

  if (defined(amrfinder_stress_reports)) {
    call tasks.build_gene_matrix as amrfinder_stress_matrix {
      input:
        tsv_reports = select_first([amrfinder_stress_reports]),
        gene_column = "Element symbol",
        source_name = "amrfinder_stress"
    }
  }

  if (defined(amrfinder_virulence_reports)) {
    call tasks.build_gene_matrix as amrfinder_virulence_matrix {
      input:
        tsv_reports = select_first([amrfinder_virulence_reports]),
        gene_column = "Element symbol",
        source_name = "amrfinder_virulence"
    }
  }

  output {
    # AbriCate
    File? abricate_matrix_yesno_tsv    = abricate_matrix.matrix_yesno_tsv
    File? abricate_matrix_binary_tsv   = abricate_matrix.matrix_binary_tsv
    File? abricate_matrix_yesno_xlsx   = abricate_matrix.matrix_yesno_xlsx
    File? abricate_matrix_binary_xlsx  = abricate_matrix.matrix_binary_xlsx
    File? abricate_concatenated_tsv    = abricate_matrix.concatenated_tsv

    # AMRFinderPlus - resistance (amr)
    File? amrfinder_amr_matrix_yesno_tsv    = amrfinder_amr_matrix.matrix_yesno_tsv
    File? amrfinder_amr_matrix_binary_tsv   = amrfinder_amr_matrix.matrix_binary_tsv
    File? amrfinder_amr_matrix_yesno_xlsx   = amrfinder_amr_matrix.matrix_yesno_xlsx
    File? amrfinder_amr_matrix_binary_xlsx  = amrfinder_amr_matrix.matrix_binary_xlsx
    File? amrfinder_amr_concatenated_tsv    = amrfinder_amr_matrix.concatenated_tsv

    # AMRFinderPlus - stress
    File? amrfinder_stress_matrix_yesno_tsv    = amrfinder_stress_matrix.matrix_yesno_tsv
    File? amrfinder_stress_matrix_binary_tsv   = amrfinder_stress_matrix.matrix_binary_tsv
    File? amrfinder_stress_matrix_yesno_xlsx   = amrfinder_stress_matrix.matrix_yesno_xlsx
    File? amrfinder_stress_matrix_binary_xlsx  = amrfinder_stress_matrix.matrix_binary_xlsx
    File? amrfinder_stress_concatenated_tsv    = amrfinder_stress_matrix.concatenated_tsv

    # AMRFinderPlus - virulence
    File? amrfinder_virulence_matrix_yesno_tsv    = amrfinder_virulence_matrix.matrix_yesno_tsv
    File? amrfinder_virulence_matrix_binary_tsv   = amrfinder_virulence_matrix.matrix_binary_tsv
    File? amrfinder_virulence_matrix_yesno_xlsx   = amrfinder_virulence_matrix.matrix_yesno_xlsx
    File? amrfinder_virulence_matrix_binary_xlsx  = amrfinder_virulence_matrix.matrix_binary_xlsx
    File? amrfinder_virulence_concatenated_tsv    = amrfinder_virulence_matrix.concatenated_tsv
  }
}
