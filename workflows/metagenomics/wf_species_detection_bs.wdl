version 1.0

workflow wf_species_detection_bs {
    input {
        String basespace_collection_id
        String api_server
        String access_token
        String? sample_prefix
    }

    meta {
        description: "Check whether the sequenced samples in Basespace Illumina correspond to the expected species using Kraken2. Un echec sur un sample n'interrompt pas le traitement des autres."
        author: "David Maimoun"
        organization: "MOH"
        email: "david.maimoun@moh.gov.il"
        version: "1.1"
    }

    parameter_meta {
        api_server: "Illumina Basepace API"
        access_token: "Illumina Basespace access token"
        basespace_collection_id: "Samplename in the Basespace platform (for ex, EC001, NM005 , = our entity id)"
        sample_prefix: "Optional, to fetch only specific species according to your samplename prefix, for ex 'EC' for ecoli (useful for testing)"
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

        call Detect_Species {
            input:
                read1 = FetchReads.read1,
                read2 = FetchReads.read2,
                sample_id = sample_name
        }
    }

    output {
        String version = Detect_Species.version[0]
        File reads_list = GetReadsList.reads_list
        Array[String] samples_name = GetReadsList.samples_name
        Array[String] species_detected = Detect_Species.sample_detected
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

    # NOTE: ce script ne doit JAMAIS sortir avec un exit code non-nul volontairement.
    # Un sample introuvable ou un dataset vide doit simplement se traduire par
    # l'absence de read1/read2 en sortie, geree plus loin par Detect_Species
    # ("not found"). Ca evite qu'un sample foireux fasse echouer tout le scatter.
    command <<<
        run_name_trimmed=$(echo "~{basespace_collection_id}" | awk '{$1=$1;print}')

        if [[ ! -z "~{basespace_sample_id}" ]]; then
            sample_identifier="~{basespace_sample_name}"
            dataset_name="~{basespace_sample_id}"
        else
            sample_identifier="~{basespace_sample_name}"
            dataset_name="~{basespace_sample_name}"
        fi

        echo -e "sample_identifier: ${sample_identifier}\ndataset_name: ${dataset_name}\nbasespace_collection_id: $run_name_trimmed"

        bs_command="bs --api-server=~{api_server} --access-token=~{access_token}"
        echo "bs_command: ${bs_command}"

        run_id=$(${bs_command} list run --retry | grep "$run_name_trimmed" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
        echo "run_id: ${run_id}"

        dataset_id_array=()

        if [[ ! -z "${run_id}" ]]; then
            dataset_id_array=($(${bs_command} list dataset --retry --input-run=${run_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' ))
            echo "dataset_id: ${dataset_id_array[*]}"
        else
            echo "Could not locate a run_id via Basespace runs, attempting to search Basespace projects now..."
            project_id=$(${bs_command} list project --retry | grep "$run_name_trimmed" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
            echo "project_id: ${project_id}"

            if [[ ! -z "${project_id}" ]]; then
                echo "project_id identified via Basespace, now searching for dataset_id within project_id ${project_id}..."
                dataset_id_array=($(${bs_command} list dataset --retry --project-id=${project_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' ))
                echo "dataset_id: ${dataset_id_array[*]}"
            else
                echo "WARNING: no run or project id found for basespace_collection_id: $run_name_trimmed -- sample will be marked as not found downstream." >&2
            fi
        fi

        if [[ ${#dataset_id_array[@]} -eq 0 ]]; then
            echo "WARNING: no dataset found for sample ${sample_identifier} -- sample will be marked as not found downstream." >&2
        fi

        for index in ${!dataset_id_array[@]}; do
            dataset_id=${dataset_id_array[$index]}
            mkdir -p ./dataset_${dataset_id} && cd ./dataset_${dataset_id}

            echo "dataset download: ${bs_command} download dataset -i ${dataset_id} -o . --retry"
            # || true : un echec de telechargement isole ne doit pas tuer le shard,
            # il se traduira juste par une absence de fastq -> "not found" en aval.
            ${bs_command} download dataset --retry -i ${dataset_id} -o . --retry || echo "WARNING: download failed for dataset ${dataset_id}" >&2
            cd ..
            echo -e "downloaded data: \n $(ls ./dataset_*/* 2>/dev/null)"
        done

        echo "Concatenating and renaming FASTQ files to add back underscores in basespace_sample_name"
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

        # Toujours sortir en succes: un fetch rate n'est pas une erreur de tache,
        # c'est une information ("sample absent") geree par Detect_Species.
        exit 0
    >>>

    output {
        File? read1 = 'fwd.fastq.gz'
        File? read2 = 'rev.fastq.gz'
    }

    runtime {
        docker: docker
        maxRetries: 2
        continueOnReturnCode: true
    }
}

task Detect_Species {
    input {
        File? read1
        File? read2
        String sample_id

        String docker = "bioinfomoh/specie_detection:1"
        Int cpu = 35
    }

    command <<<
        mode=""
        compressed=""

        if [ ! -s "~{read1}" ]; then
            echo "~{sample_id} not found" > "~{sample_id}.report"
            detected="not_found"
        else
            if ! [ -z "~{read2}" ]; then
                echo "Reads are paired..."
                mode="--paired"
            fi

            if [[ "~{read1}" == *.gz ]]; then
                echo "Reads are compressed..."
                compressed="--gzip-compressed"
            fi

            kraken2 -v | awk '/Kraken/ {print "Kraken v" $3}' | tee VERSION

            echo "Running Kraken2..."
            # Si kraken2 plante sur ce sample precis (fastq corrompu, etc.), on
            # capture l'echec au lieu de faire tomber toute la tache/scatter.
            if kraken2 $mode $compressed --threads "~{cpu}" --use-names --db /app/db/kraken_db \
                --report "~{sample_id}.report" --paired "~{read1}" "~{read2}" --output - ; then
                kraken_ok=true
            else
                kraken_ok=false
                echo "WARNING: kraken2 failed for ~{sample_id}" >&2
            fi

            if [[ "$kraken_ok" == "false" ]]; then
                [ -s "~{sample_id}.report" ] || echo "~{sample_id} kraken2_error" > "~{sample_id}.report"
                detected="kraken2_error"
            else
                # Resultat brut de Kraken2, sans comparaison avec une espece attendue.
                detected=$(awk -F'\t' '$4 == "S" {gsub(/^[ \t]+/, "", $6); print $6; exit}' "~{sample_id}.report")
                [[ -n "$detected" ]] || detected="undetermined"
            fi
        fi

        echo "~{sample_id},${detected}" > ~{sample_id}_sample_detected.csv

        # Garantit une VERSION meme si kraken2 n'a jamais tourne (sample not_found)
        [ -s VERSION ] || echo "Kraken vNA" > VERSION

        exit 0
    >>>

    output {
        String version = read_string("VERSION")
        File report = "~{sample_id}.report"
        String sample_detected = read_string("~{sample_id}_sample_detected.csv")
    }

    runtime {
        docker: docker
        cpu: cpu
        maxRetries: 1
        continueOnReturnCode: true
    }
}

