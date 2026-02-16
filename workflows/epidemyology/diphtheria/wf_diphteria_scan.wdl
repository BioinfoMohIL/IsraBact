version 1.0

import "../../../tasks/epidemiology/diphtheria/task_diphtoscan.wdl" as task_diphtoscan
import "../../../tasks/task_fail.wdl" as task_fail

workflow wf_diphtheria_scan {
    meta {
        description: "diphtOscan is a tool to search genomic assemblies of Corynebacterium diphtheriae and other species of the Corynebacterium diphtheriae species complex (CdSC) - (https://gitlab.pasteur.fr/BEBP/diphtoscan) ."
        author: "David Maimoun"
    }

    input {
        Array[File]? assemblies
        File? assembly

        Boolean mlst = true
        Boolean tox = true
        Boolean res_vir = true
        Boolean integron = true
        Boolean extend_genotyping = false
        Boolean tree = true
        Boolean update_db = true

        Int threads = 16
    }

    Array[File] input_files = flatten([
        select_first([assemblies, []]),
        if defined(assembly) then [select_first([assembly])] else []
    ])

    if (length(input_files) == 0) {
        call task_fail.fail { 
            input: message = "Must provide assembliy or assemblies" 
        }
    }

    call task_diphtoscan.diphtoscan {
        input:
            my_input = input_files,
            mlst = mlst,
            tox = tox,
            res_vir = res_vir,
            extend_genotyping = extend_genotyping,
            integron = integron,
            tree = tree,
            threads = threads,
            update_db = update_db
    }


    output {
        File diphtoscan_results = select_first([diphtoscan.results, fail.fail_logs])
        File? diphtoscan_tree_results = diphtoscan.tree_results

    }
}


