version 1.0

import "../../../tasks/epidemiology/diphtheria/task_diphtheria_scan.wdl" as task_diphtheria_scan

workflow wf_diphtheria_scan {
    input {
        Array[File] assemblies
        Boolean mlst = true
        Boolean tox = true
        Boolean res_vir = true
        Boolean extend_genotyping = false
        Boolean tree = true
        Boolean integron = false
        Int threads = 4
    }

    call task_diphtheria_scan.run_diphtOscan {
        input:
            assemblies = assemblies,
            mlst = mlst,
            tox = tox,
            res_vir = res_vir,
            extend_genotyping = extend_genotyping,
            integron = integron,
            tree = tree,
            threads = threads
    }

    output {
        Array[File] results = run_diphtOscan.results
    }
}


