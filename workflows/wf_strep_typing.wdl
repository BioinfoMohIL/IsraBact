version 1.0

import "../tasks/tasks_versioning.wdl" as versioning
import "../tasks/typing/streptococcus/task_seroba_v2" as seroba_v2

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

    call seroba_v2 as seroba {
        input:
            read1 = read1,
            read2 = read2,
            samplename = samplename
    }

    call versioning.version_capture {
        input:
    }

    output {
        String basespace_fetch_version = version_capture.phb_version
        String basespace_fetch_analysis_date = version_capture.date

        String seroba_v2_version = seroba.seroba_version
        String seroba_v2_docker = seroba.seroba_docker
        String seroba_v2_serotype = seroba.seroba_serotype
        String? seroba_v2_ariba_serotype = seroba.seroba_ariba_serotype
        String? seroba_v2_ariba_identity = seroba.seroba_ariba_identity
        File? seroba_v2_details = seroba.seroba_details
    }

}