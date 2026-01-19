version 1.0

import "tasks/utilities/task_matrixer.wdl" as task_matrixer

workflow wf_matrixer {
    input {
        File input_tsv
        String? genes
    }

    call task_matrixer.build_matrix {
        input:
            amr_tsv = input_tsv,
            genes = genes
    }

    output {
        File matrixer_matrix = build_matrix.output_matrix
    }
}




