version 1.0

import "../../tasks/data_visualization/task_build_gene_matrix.wdl" as tasks

## 4 optional inputs, all independent from one another:
##   - abricate_results_tsvs             : AbriCate reports (virulence genes)
##   - amrfinder_amr_reports        : AMRFinderPlus reports - resistance genes (amr_report)
##   - amrfinder_stress_reports     : AMRFinderPlus reports - stress genes (stress_report)
##   - amrfinder_virulence_reports  : AMRFinderPlus reports - virulence genes (virulence_report)
##
## Only the inputs that are provided get processed: if an input is omitted, the
## corresponding task is not run and its outputs stay empty. This makes it possible,
## for example, to generate only a virulence matrix without having to supply the
## other 3 inputs.
##
## For each input provided, produces:
##   - YES/NO matrix   (TSV + XLSX)
##   - binary matrix   (TSV + XLSX)
##   - a file concatenating the raw input reports (with a Key column)


workflow wf_amr_gene_matrix {
  input {
    Array[File]? abricate_results_tsvs
    Array[File]? amrfinder_amr_reports
    Array[File]? amrfinder_stress_reports
    Array[File]? amrfinder_virulence_reports
  }

  if (defined(abricate_results_tsvs)) {
    call tasks.build_gene_matrix as abricate_matrix {
      input:
        tsv_reports = select_first([abricate_results_tsvs]),
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
    # File? abricate_matrix_yesno_tsv    = abricate_matrix.matrix_yesno_tsv
    # File? abricate_matrix_binary_tsv   = abricate_matrix.matrix_binary_tsv
    File? abricate_matrix_yesno_xlsx   = abricate_matrix.matrix_yesno_xlsx
    File? abricate_matrix_binary_xlsx  = abricate_matrix.matrix_binary_xlsx
    File? abricate_concatenated_tsv    = abricate_matrix.concatenated_tsv

    # AMRFinderPlus - resistance (amr)
    # File? amrfinder_amr_matrix_yesno_tsv    = amrfinder_amr_matrix.matrix_yesno_tsv
    # File? amrfinder_amr_matrix_binary_tsv   = amrfinder_amr_matrix.matrix_binary_tsv
    File? amrfinder_amr_matrix_yesno_xlsx   = amrfinder_amr_matrix.matrix_yesno_xlsx
    File? amrfinder_amr_matrix_binary_xlsx  = amrfinder_amr_matrix.matrix_binary_xlsx
    File? amrfinder_amr_concatenated_tsv    = amrfinder_amr_matrix.concatenated_tsv

    # AMRFinderPlus - stress
    # File? amrfinder_stress_matrix_yesno_tsv    = amrfinder_stress_matrix.matrix_yesno_tsv
    # File? amrfinder_stress_matrix_binary_tsv   = amrfinder_stress_matrix.matrix_binary_tsv
    File? amrfinder_stress_matrix_yesno_xlsx   = amrfinder_stress_matrix.matrix_yesno_xlsx
    File? amrfinder_stress_matrix_binary_xlsx  = amrfinder_stress_matrix.matrix_binary_xlsx
    File? amrfinder_stress_concatenated_tsv    = amrfinder_stress_matrix.concatenated_tsv

    # AMRFinderPlus - virulence
    # File? amrfinder_virulence_matrix_yesno_tsv    = amrfinder_virulence_matrix.matrix_yesno_tsv
    # File? amrfinder_virulence_matrix_binary_tsv   = amrfinder_virulence_matrix.matrix_binary_tsv
    File? amrfinder_virulence_matrix_yesno_xlsx   = amrfinder_virulence_matrix.matrix_yesno_xlsx
    File? amrfinder_virulence_matrix_binary_xlsx  = amrfinder_virulence_matrix.matrix_binary_xlsx
    File? amrfinder_virulence_concatenated_tsv    = amrfinder_virulence_matrix.concatenated_tsv
  }
}
