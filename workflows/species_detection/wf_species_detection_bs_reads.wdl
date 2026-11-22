version 1.0

import "wf_species_detection.wdl" as detection

# BaseSpace version of the species detection: it only fetches the reads from BaseSpace, then runs
# wf_species_detection (same detection, same match logic, same report/summary/plate check as in Terra).

workflow wf_species_detection_bs {
    input {
        String basespace_collection_id
        String api_server
        String access_token
        String? sample_prefix

        File?   sending_file          # if set, the plate check xlsx is also produced
        Boolean use_ml = false        # false = kraken2, true = ML model
        File?   ml_model
        Float?  ml_threshold
        String  docker = "ghcr.io/bioinfomohil/species_detection:2.0"
    }

    meta {
        description: "Check whether the sequenced samples in Basespace Illumina correspond to the expected species (kraken2 or ML)."
        author: "David Maimoun"
        organization : "MOH"
        email: "david.maimoun@moh.gov.il"
        version: "2.0"
    }

    parameter_meta {
        api_server: "Illumina Basepace API"
        access_token: "Illumina Basespace access token"
        basespace_collection_id: "Samplename in the Basespace platform (for ex, EC001, NM005 , = our entity id)"
        sample_prefix: "Optional, to fetch only specific species according to your samplename prefix, for ex 'EC' for ecoli (useful for testing)"
        sending_file: "Optional, the lab sending file (.xlsx). If set, a plate check xlsx is also produced (one sheet per plate). A failure of the plate check does not fail the workflow."
        use_ml: "false (default) = kraken2, true = ML model"
    }

    call GetReadsList {
        input:
            basespace_collection_id = basespace_collection_id,
            access_token = access_token,
            api_server = api_server,
            sample_prefix = select_first([sample_prefix, ""])
    }

    scatter (sample_name in GetReadsList.samples_name) {
        call FetchReads {
          input:
            basespace_sample_name = sample_name,
            basespace_collection_id = basespace_collection_id,
            api_server = api_server,
            access_token = access_token
        }

        # only the samples for which both reads were downloaded go to the detection
        if (defined(FetchReads.read1) && defined(FetchReads.read2)) {
            String found_id = sample_name
            File   found_r1 = select_first([FetchReads.read1])
            File   found_r2 = select_first([FetchReads.read2])
        }
        if (!(defined(FetchReads.read1) && defined(FetchReads.read2))) {
            String missing_id = sample_name
        }
    }

    call detection.wf_species_detection as detect {
        input:
            read1 = select_all(found_r1),
            read2 = select_all(found_r2),
            sample_ids = select_all(found_id),
            sending_file = sending_file,
            use_ml = use_ml,
            ml_model = ml_model,
            ml_threshold = ml_threshold,
            docker = docker
    }

    output {
        File          reads_list               = GetReadsList.reads_list
        Array[String] samples_name             = GetReadsList.samples_name
        Array[String] samples_without_reads    = select_all(missing_id)   # not in the report
        String        sp_detection_scanner     = detect.sp_detection_scanner
        File          sp_detection_report      = detect.sp_detection_report
        String        sp_detection_summary     = detect.sp_detection_summary
        File?         sp_detection_plate_check = detect.sp_detection_plate_check
        String?       sp_detection_plate_status = detect.sp_detection_plate_status
    }
}

task GetReadsList {
    input {  
        String basespace_collection_id
        String api_server 
        String access_token
        String? sample_prefix
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"

    }

    command <<<     
        bs project content --name ~{basespace_collection_id} \
            --api-server=~{api_server} \
            --access-token=~{access_token} \
            --retry > reads_list.txt

        # Create a a samples name list from the reads_list
        # -> Fetch the name form the R1 read, and remove if got a 'Undetermined' read
        if [ -z ~{sample_prefix} ]; then
            grep -o "[A-Za-z0-9_-]*_S[0-9]*_L[0-9]*_R1_[0-9]*\.fastq\.gz" reads_list.txt \
            | sed 's/_S[0-9]*_L[0-9]*_R1_.*\.fastq\.gz//' \
            | grep -v "^Undetermined$" \
            > samples_name.txt
        else
            grep -o "~{sample_prefix}[A-Za-z0-9_-]*_S[0-9]*_L[0-9]*_R1_[0-9]*\.fastq\.gz" reads_list.txt \
            | sed 's/_S[0-9]*_L[0-9]*_R1_.*\.fastq\.gz//' \
            > samples_name.txt
        fi
    >>>

    output {
        File reads_list = "reads_list.txt"
        Array[String] samples_name = read_lines("samples_name.txt")
    }

    runtime {
        docker: docker
        preemptible: 1
  }
}

task FetchReads {
    input {
        String basespace_sample_name
        String? basespace_sample_id   
        String basespace_collection_id
        String api_server 
        String access_token
     
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"
    }

    command <<<
        # Trim basespace_collection_id (run_name) because spaces cause problem
        run_name_trimmed=$(echo "~{basespace_collection_id}" | awk '{$1=$1;print}')

        # set basespace name and id variables
        if [[ ! -z "~{basespace_sample_id}" ]]; then
            sample_identifier="~{basespace_sample_name}"
            dataset_name="~{basespace_sample_id}"
        else
            sample_identifier="~{basespace_sample_name}"
            dataset_name="~{basespace_sample_name}"
        fi
    
        # print all relevant input variables to stdout
        echo -e "sample_identifier: ${sample_identifier}\ndataset_name: ${dataset_name}\nbasespace_collection_id: $run_name_trimmed"
        
        #Set BaseSpace comand prefix
        bs_command="bs --api-server=~{api_server} --access-token=~{access_token}"
        echo "bs_command: ${bs_command}"

        #Grab BaseSpace Run_ID from given BaseSpace Run Name
        run_id=$(${bs_command} list run --retry | grep "$run_name_trimmed" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
        echo "run_id: ${run_id}" 
        
        if [[ ! -z "${run_id}" ]]; then 
            #Grab BaseSpace Dataset ID from dataset lists within given run 
            dataset_id_array=($(${bs_command} list dataset --retry --input-run=${run_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' )) 
            echo "dataset_id: ${dataset_id_array[*]}"
        
        else 
            #Try Grabbing BaseSpace Dataset ID from project name
            echo "Could not locate a run_id via Basespace runs, attempting to search Basespace projects now..."
            
            project_id=$(${bs_command} list project --retry | grep "$run_name_trimmed" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
            
            echo "project_id: ${project_id}" 

            if [[ ! -z "${project_id}" ]]; then 
                echo "project_id identified via Basespace, now searching for dataset_id within project_id ${project_id}..."
                dataset_id_array=($(${bs_command} list dataset --retry --project-id=${project_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' ))  
                echo "dataset_id: ${dataset_id_array[*]}"
            else       
                echo "No run or project id found associated with input basespace_collection_id: $run_name_trimmed" >&2
                exit 1
            fi      
        fi

        #Download reads by dataset ID
        for index in ${!dataset_id_array[@]}; do
            dataset_id=${dataset_id_array[$index]}
            mkdir ./dataset_${dataset_id} && cd ./dataset_${dataset_id}
        
            echo "dataset download: ${bs_command} download dataset -i ${dataset_id} -o . --retry"
            ${bs_command} download dataset --retry -i ${dataset_id} -o . --retry && cd ..
            echo -e "downloaded data: \n $(ls ./dataset_*/*)"
        done

        # rename FASTQ files to add back in underscores that Illumina/Basespace changed into hyphens
        echo "Concatenating and renaming FASTQ files to add back underscores in basespace_sample_name"
        # setting a new bash variable to use for renaming during concatenation of FASTQs
        # SAMPLENAME_HYPHEN_INSTEAD_OF_UNDERSCORES=$(echo $sample_identifier | sed 's|_|-|g' | sed 's|\.|-|g')

        # echo $SAMPLENAME_HYPHEN_INSTEAD_OF_UNDERSCORES > 'sample_id.txt'

        echo $sample_identifier > 'sample_id.txt'
 

        for fwd_read in ./dataset_*/${sample_identifier}_*R1_*.fastq.gz; do
            if [[ -s $fwd_read ]]; then
                read1_name=$(basename "$fwd_read")

                echo ${read1_name} > read1_name.txt
                cat $fwd_read      > fwd.fastq.gz
                
            fi
        done

        for rev_read in ./dataset_*/${sample_identifier}_*R2_*.fastq.gz; do
            if [[ -s $rev_read ]]; then
                read2_name=$(basename "$rev_read")

                echo ${read2_name} > read2_name.txt
                cat $rev_read      > rev.fastq.gz
            fi
        done
    >>>

    output {
        File? read1          = 'fwd.fastq.gz'
        File? read2          = 'rev.fastq.gz'
    }

    runtime {
        docker: docker
        maxRetries: 1
        
  }
}
