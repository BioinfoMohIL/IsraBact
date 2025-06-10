version 1.0

import "../tasks/tasks_versioning.wdl" as versioning_task
import '../tasks/typing/legionella/task_elgato.wdl' as elgato

workflow sbt_analysis {
    meta {
        description: "SBT Legionella pneumophila isolates using the elgato tool."
        author: "David Maimoun"
        organization: "MOH Jerusalem"
    }

    input {
        String samplename
        File read1 
        File read2
        String docker = "staphb/elgato:1.21.2"
    }

    
    call elgato.elgato_reads as elgato_reads{
        input:
        read1       = read1,
        read2       = read2,
        samplename  = samplename,
        docker      = docker
    }

    call version_capture {
        input:
    }

    output {
        String basespace_fetch_version = version_capture.phb_version
        String basespace_fetch_analysis_date = version_capture.date

        String sbt_elgato_version = elgato_reads.elgato_version
        String sbt = elgato_reads.sbt
        File sbt_possible_sts = elgato_reads.possible_mlsts
        File sbt_inter_out = elgato_reads.intermediate_outputs
        File sbt_alleles = elgato_reads.alleles
    }
}



