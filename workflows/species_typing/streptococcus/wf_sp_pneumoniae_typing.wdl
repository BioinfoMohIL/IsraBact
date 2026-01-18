version 1.0

import "../../../tasks/species_typing/streptococcus/task_seroba_v2.wdl" as seroba

workflow strep_typing {
    meta {
        description: "Streptococcus Serotying using Seroba v2.0.4."
        author: "David Maimoun"
        organization: "MOH Jerusalem"
    }

    input {
        File read1
        File read2 
        String samplename 
    }

    call seroba.seroba_v2  {
        input:
            read1 = read1,
            read2 = read2,
            samplename = samplename
    }

    output {
        String seroba_v2_version = seroba_v2.seroba_version
        String seroba_v2_docker = seroba_v2.seroba_docker
        String seroba_v2_serotype = seroba_v2.seroba_serotype
        String? seroba_v2_ariba_serotype = seroba_v2.seroba_ariba_serotype
        String? seroba_v2_ariba_identity = seroba_v2.seroba_ariba_identity
        File? seroba_v2_details = seroba_v2.seroba_details
    }

}