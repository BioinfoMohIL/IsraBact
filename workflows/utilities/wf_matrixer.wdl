version 1.0

import "../../tasks/utilities/task_matrixer.wdl" as task_matrixer

workflow wf_matrixer {
    input {
        File input_tsv
        String scope
    }

    call task_matrixer.build_matrix {
        input:
            input_tsv = input_tsv,
            scope = scope
    }

    output {
        File matrixer_matrix = build_matrix.output_matrix
        String matrixer_scope = scope
    }
}




