version 1.0

import "../../tasks/data_visualization/task_build_gene_matrix.wdl" as tasks

## wf_amr_gene_matrix.wdl
##
## Prend en entree, pour un lot de souches :
##   - abricate_reports  : Array[File] des rapports AbriCate (genes de virulence)
##   - amrfinder_reports : Array[File] des rapports AMRFinderPlus (resistance/stress/virulence)
##
## Pour chaque source, produit 2 matrices "souche x gene" (TSV + XLSX) :
##   - YES / NO  (presence/absence)
##   - 0 / 1     (binaire)

workflow wf_amr_gene_matrix {
  input {
    Array[File] abricate_results
    Array[File] amrfinder_reports
  }

  call tasks.build_gene_matrix as abricate_matrix {
    input:
      tsv_reports = abricate_results,
      gene_column = "GENE",
      source_name = "abricate"
  }

  call tasks.build_gene_matrix as amrfinder_matrix {
    input:
      tsv_reports = amrfinder_reports,
      gene_column = "Element symbol",
      source_name = "amrfinder"
  }

  output {
    File amr_gene_matrix_abricate_yesno_tsv    = abricate_matrix.matrix_yesno_tsv
    File amr_gene_matrix_abricate_binary_tsv   = abricate_matrix.matrix_binary_tsv
    File amr_gene_matrix_abricate_yesno_xl     = abricate_matrix.matrix_yesno_xlsx
    File amr_gene_matrix_abricate_binary_xl    = abricate_matrix.matrix_binary_xlsx

    File amr_gene_matrix_amrfinder_yesno_tsv   = amrfinder_matrix.matrix_yesno_tsv
    File amr_gene_matrix_amrfinder_binary_tsv  = amrfinder_matrix.matrix_binary_tsv
    File amr_gene_matrix_amrfinder_yesno_xl    = amrfinder_matrix.matrix_yesno_xlsx
    File amr_gene_matrix_amrfinder_binary_xl   = amrfinder_matrix.matrix_binary_xlsx
  }
}
