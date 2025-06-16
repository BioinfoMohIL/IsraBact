version 1.0

import "../tasks/utilities/data_import/task_basespace_cli.wdl" as basespace
import "../tasks/metagenomics/task_detect_species.wdl" as metagenomics

workflow species_detection {
    meta {
        description: "Detect your species using Kraken2."
        author: "David Maimoun (The Codon Bleu)"
        email: "thecodonbleu@outlook.com"
    }

    input {
        String basespace_collection_id 
        String api_server
        String access_token
        String? sample_prefix
    }

    call basespace.get_reads_list as get_reads_list {
        input:
            basespace_collection_id = basespace_collection_id,
            access_token = access_token,
            api_server = api_server,
            sample_prefix = select_first([sample_prefix, ""])

    }

    scatter(sample_name in get_reads_list.samples_name) {
        call basespace.fetch_bs as fetch_bs {
          input:
            sample_name = sample_name,
            basespace_sample_name = sample_name,
            basespace_collection_id = basespace_collection_id,
            api_server = api_server,
            access_token = access_token
        }

        call metagenomics.detect_species  {
            input:
                read1 = fetch_bs.read1,
                read2 = fetch_bs.read2,
                sample_id = sample_name
        }
        
    }

    call merge_reports {
        input:
            species_detected_list = detect_species.sample_detected
    }

  
    output {
    
        File reads_list = get_reads_list.reads_list
        Array[String] samples_name = get_reads_list.samples_name
        File species_detected = merge_reports.species_detected
    
    }

}
  



task merge_reports {
    input {
        Array[String] species_detected_list
    }

    command <<<
        echo "Sample,Detected,Match" > species_detected.csv
        echo "~{sep='\n' species_detected_list}" >> species_detected.csv
        
    >>>

    output {
        File species_detected = "species_detected.csv"
    }

    runtime {
        docker: "ubuntu:20.04"

    }
}