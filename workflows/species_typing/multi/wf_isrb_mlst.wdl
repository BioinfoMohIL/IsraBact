version 1.0

import "../../../tasks/species_typing/multi/task_isrb_mlst.wdl" as task_mlst

workflow wf_mlst {
    input {
        File assembly_fasta
        String taxon
    }

    call task_mlst.mlst as mlst {
        input:
            assembly_fasta = assembly_fasta,
            taxon          = taxon
    }

    output {
        File   isrb_mlst_result   = mlst.result
        String isrb_mlst_db       = mlst.db_used
        String isrb_mlst_alleles  = mlst.alleles
        String isrb_mlst_st       = mlst.st
    }
}
