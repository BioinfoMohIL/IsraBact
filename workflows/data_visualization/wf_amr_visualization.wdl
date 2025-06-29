version 1.0

workflow run_amr_heatmap {

  input {
    File amr_input
    Boolean? remove_null = true
  }

  call amr_heatmap {
    input:
      amr_input = amr_input,
      remove_null = remove_null,
  }

  output {
    File plots = amr_heatmap.plots
  }
}

task amr_heatmap {

  input {
    File amr_input           # CSV or XLSX input file
    Boolean? remove_null       # Optional flag to remove NA rows
  }

  command <<<
    mkdir plots
    amr_analysis \
        --amr_input ~{amr_input} \
        --output plots \
        ~{true='--remove_null' false='' remove_null}

    zip -r plots.zip plots

  >>>

  output {
    File plots = "plots.zip"

  }

  runtime {
    docker: "bioinfomoh/data_analysis_r_tools:1"
  }
}
