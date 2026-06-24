version 1.0

task fetch_species_schema_adapted {
    input {
        String species
        Int disk_gb = 50
        Int memory_gb = 4

        # Optional custom archive provided by user
        File? schema_adapted_zip
    }

    command <<<
        set -euo pipefail

        SPECIES="~{species}"

        echo "[info] Requested species: $SPECIES"

        if [ "~{defined(schema_adapted_zip)}" = "true" ]; then
            echo "[info] Using provided schema archive"

            cp "~{schema_adapted_zip}" schema_adapted_archive

        else
            echo "[info] No custom archive provided, using default GCS schema"

            case "$SPECIES" in
                Campylobacter)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/campylobacter_jejuni/Campylobacter_jejuni_wgMLST_2025-07-31T16_53_57.662910.zip"
                    ;;
                Ecoli)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/ecoli/Ecoli_wgMLST.zip"
                    ;;
                Salmonella)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/salmonella/Salmonella_wgMLST.zip"
                    ;;
                Streptococcus)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/streptococcus/Streptococcus_wgMLST.zip"
                    ;;
                Listeria)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/listeria/Listeria_wgMLST.zip"
                    ;;
                *)
                    echo "[error] Unsupported species: $SPECIES"
                    exit 1
                    ;;
            esac

            echo "[info] Downloading archive from:"
            echo "$GCS_PATH"

            cp "$GCS_PATH" schema_adapted_archive.zip

            mv schema_adapted_archive.zip schema_adapted_archive
        fi


        mkdir -p schema_adapted

        unzip -q schema_adapted_archive -d schema_adapted

    >>>

    output {
        Array[File] schema_adapted_files = glob("schema_adapted/**/*")
    }

    runtime {
        docker: "python:3.11"
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}

task allele_calling {
    input {
        String specie_prefix
        Array[File] assemblies
        File? schema_adapted_zip
        File? assemblies_zipped
        Int     allele_call_mode               
        Float   allele_call_bsr              
        Int     allele_call_minimum_length    
        Int     allele_call_translation_table 
        Float   allele_call_size_threshold    
        Boolean allele_call_cds_input        
        Boolean allele_call_no_inferred       
        
        Int cpu = 30
    }

    command <<<
        set -e

        PREFIX="~{specie_prefix}"
        PREFIX=$(echo "$PREFIX" | tr '[:lower:]' '[:upper:]')
        
        echo "[info] Requested species: $SPECIES"

        if [ "~{defined(schema_adapted_zip)}" = "true" ]; then
            echo "[info] Using provided schema archive"
            cp "~{schema_adapted_zip}" schema_adapted_archive

        else
            case "$PREFIX" in
                CA)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/campylobacter_jejuni/Campylobacter_jejuni_wgMLST_2025-07-31T16_53_57.662910.zip"
                    ;;
                EC|F-EC)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/escherichia_coli/Escherichia_coli_wgMLST_2025-08-01T13_27_15.392572.zip"
                    ;;
                SA|SO)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/salmonella_enterica/Salmonella_enterica_wgMLST_2025-08-01T20_14_39.891230.zip"
                    ;;
                ST)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/streptococcus_pyogenes/Streptococcus_pyogenes_wgMLST_2025-07-23T16_14_20.901319.zip"
                    ;;
                NM|M)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/neisseria_meningitidis/Neisseria_meningitidis_cgMLST_2025-08-01T16_08_32.955509.zip"
                    ;;
                 LG)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/legionella_pneumophila/Legionella_pneumophila_wgMLST_2026-02-18T15_35_04.061183.zip"
                    ;;
                LF|LC)
                    GCS_PATH="gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/chewbbaca/listeria/Listeria_wgMLST.zip"
                    ;;
                *)
                    echo "[error] Unsupported species prefix: $PREFIX"
                    exit 1
                    ;;

            esac

            echo "[info] Downloading archive from:"
            echo "$GCS_PATH"

            gsutil cp "$GCS_PATH" schema_adapted_archive.zip

            mv schema_adapted_archive.zip schema_adapted_archive
           
        fi

        mkdir -p schema_adapted
        unzip -q schema_adapted_archive -d schema_adapted

        mkdir -p assemblies_files

        echo "[2] Preparing assemblies input"
        # Cas 1 : assemblies_zipped 
        if [ -s "~{assemblies_zipped}" ]; then
            echo "Decompressing assemblies_zipped"
            tmp_dir=$(mktemp -d)

            if [[ "$(basename "~{assemblies_zipped}")" == *.zip ]]; then
                unzip -q "~{assemblies_zipped}" -d "$tmp_dir"
            else
                tar -xf "~{assemblies_zipped}" -C "$tmp_dir"
            fi

            echo "Moving files to assemblies_files"
            find "$tmp_dir" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) -exec mv {} assemblies_files/ \;
            rm -rf "$tmp_dir"
        fi
        
        # Cas 2 : assemblies array files
        if [ "~{sep=' ' assemblies}" != "" ]; then
            echo "Copying assemblies list"
            for f in ~{sep=' ' assemblies}; do
                cp "$f" assemblies_files/
            done
        fi

        echo "[3] Running AlleleCall (chewBBACA v3 - Fast Mode)"
        chewBBACA.py AlleleCall \
                    -i assemblies_files \
                    -g schema_adapted \
                    -o calling \
                    --cpu ~{cpu} \

                    # --mode ~{allele_call_mode} \
                    # --bsr ~{allele_call_bsr} \
                    # --l ~{allele_call_minimum_length} \
                    # --t ~{allele_call_translation_table} \
                    # --st ~{allele_call_size_threshold} \
                    # ~{true='--cds' false='' allele_call_cds_input} \
                    # ~{true='--no_inferred' false='' allele_call_no_inferred} \

        echo "[4] Cleaning sample names"
        results_file=$(find calling -name "results_alleles.tsv" | head -n 1)
        
        awk -F'\t' 'BEGIN{OFS="\t"} {
            if (NR == 1) {
                print
            } else {
                name = $1
                sub(/\.[^.]+$/, "", name)
                split(name, parts, /[_-]/)
                $1 = parts[1]
                print
            }
        }' "$results_file" > alleles_cleaned.tsv

        echo "[5] Generating V3 Interactive Evaluation Report"
        chewBBACA.py AlleleCallEvaluator \
                    -i calling \
                    -g schema_adapted \
                    -o evaluation_report

        echo "Compressing outputs for Terra..."
        tar -czf schema_adapted.tar.gz schema_adapted
        tar -czf calling.tar.gz calling
        tar -czf evaluation_report.tar.gz evaluation_report
    >>>

    output {
        File calling_dir = "calling.tar.gz"
        File calling_alleles_cleaned = "alleles_cleaned.tsv"
        File calling_evaluation_report = "evaluation_report.tar.gz"
    }

    runtime {
        docker: "bioinfomoh/chewbbaca:3.3.10"
        cpu : "~{cpu}"
        memory: "32 GB"
        disks: "local-disk 100 HDD"
    }
}

task extract_cgmlst {
    input {
        File cleaned_results
        Array[Float] cgmlst_thresholds
    }

    command <<<
        echo "[1] Extraction of core genome loci"
        chewBBACA.py ExtractCgMLST -i ~{cleaned_results} -o visualization --t ~{sep=' ' cgmlst_thresholds}

        tar -czf visualization.tar.gz visualization
    >>>

    output {
        File visualization_zip = "visualization.tar.gz"
        Array[File] matrixes = glob("visualization/cgMLST*.tsv")

    }

    runtime {
        docker: "ummidock/chewbbaca:v3.3.10"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 20 HDD"
    }
}

task prep_external_schema {
    input {
        File    external_schema_zip     # dossier de .fasta externes zippé
        File?   training_file           # .trn prodigal (recommandé)
        File?   genes_list              # sous-ensemble de loci à adapter

        Float   bsr               = 0.6
        Int     minimum_length    = 0
        Int     translation_table = 11
        Float   size_threshold    = 0.2
        Boolean size_filter       = false
        Int     cpu               = 8

        Int     disk_gb   = 100
        Int     memory_gb = 16
    }

    command <<<
        set -euo pipefail

        echo "[decompress] external_schema_zip"
        mkdir -p external_schema
        case "~{external_schema_zip}" in
            *.zip)      unzip  "~{external_schema_zip}" -d external_schema ;;
            *.tar.gz|*.tgz) tar -xzf "~{external_schema_zip}" -C external_schema ;;
            *.tar)      tar -xf  "~{external_schema_zip}" -C external_schema ;;
            *) echo "Format non supporté"; exit 1 ;;
        esac

        # Si le zip contient un sous-dossier unique, on descend dedans
        SCHEMA_DIR=$(find external_schema -mindepth 1 -maxdepth 1 -type d | head -1)
        [ -z "$SCHEMA_DIR" ] && SCHEMA_DIR="external_schema"

        PTF_ARG=""
        [ -s "~{training_file}" ] && PTF_ARG="--ptf ~{training_file}"

        GL_ARG=""
        [ -s "~{genes_list}" ]    && GL_ARG="--gl ~{genes_list}"

        SF_ARG=""
        [ "~{size_filter}" = "true" ] && SF_ARG="--size-filter"

        echo "[run] PrepExternalSchema"
        chewBBACA.py PrepExternalSchema \
            -g "$SCHEMA_DIR" \
            -o schema_adapted \
            --bsr ~{bsr} \
            --l ~{minimum_length} \
            --t ~{translation_table} \
            --st ~{size_threshold} \
            --cpu ~{cpu} \
            $PTF_ARG \
            $GL_ARG \
            $SF_ARG

        echo "[compress] schema_adapted → schema_adapted.tar.gz"
        tar -czf schema_adapted.tar.gz schema_adapted
    >>>

    output {
        File schema_adapted_tar = "schema_adapted.tar.gz"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}

task remove_genes {
    input {
        File    results_alleles     # TSV profiles
        File    genes_list          # list of loci to remove (one per line)
        Boolean inverse = false     # if true: keep listed genes, remove others

        Int     disk_gb   = 20
        Int     memory_gb = 4
    }

    command <<<
        set -euo pipefail

        INV_ARG=""; [ "~{inverse}" = "true" ] && INV_ARG="--inverse"

        echo "[run] RemoveGenes"
        chewBBACA.py RemoveGenes \
            -i "~{results_alleles}" \
            -g "~{genes_list}" \
            -o profiles_filtered.tsv \
            $INV_ARG
    >>>

    output {
        File profiles_filtered = "profiles_filtered.tsv"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        1
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 2
    }
}

task get_alleles {
    input {
        File    results_alleles
        String  schema_dir            # GCP path
        File?   genes_list
        Boolean distinct          = false
        Boolean translate         = false
        Int     translation_table = 11
        Int     cpu               = 4

        Int     disk_gb   = 50
        Int     memory_gb = 8
    }

    command <<<
        set -euo pipefail

        echo "[gsutil] Download schema"
        mkdir -p schema_local
        gsutil -m cp -r "~{schema_dir}*" schema_local/

        GL_ARG="";    [ -s "~{genes_list}" ] && GL_ARG="--gl ~{genes_list}"
        DIST_ARG="";  [ "~{distinct}"   = "true" ] && DIST_ARG="--distinct"
        TR_ARG="";    [ "~{translate}"  = "true" ] && TR_ARG="--translate"

        echo "[run] GetAlleles"
        chewBBACA.py GetAlleles \
            -i "~{results_alleles}" \
            -g schema_local \
            -o alleles_output \
            --cpu ~{cpu} \
            --ta ~{translation_table} \
            $GL_ARG \
            $DIST_ARG \
            $TR_ARG

        tar -czf alleles_output.tar.gz alleles_output
    >>>

    output {
        File alleles_tar = "alleles_output.tar.gz"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}

task join_profiles {
    input {
        Array[File] profile_files   # list of results_alleles.tsv to join

        Int     disk_gb   = 20
        Int     memory_gb = 4
    }

    command <<<
        set -euo pipefail

        # Écrire la liste des fichiers dans un fichier texte
        echo "[run] JoinProfiles"
        chewBBACA.py JoinProfiles \
            -p ~{sep=' ' profile_files} \
            -o joined_profiles.tsv
    >>>

    output {
        File joined_profiles = "joined_profiles.tsv"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        1
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 2
    }
}

task uniprot_finder {
    input {
        String  schema_dir             # GCP path
        File?   protein_table          # cds_coordinates.tsv from CreateSchema
        File?   genes_list
        String? taxa                   # ex: "Escherichia coli"
        Float   bsr             = 0.6
        Int     proteome_matches= 1
        Boolean no_sparql       = false
        Int     cpu             = 4

        Int     disk_gb   = 50
        Int     memory_gb = 8
    }

    command <<<
        set -euo pipefail

        echo "[gsutil] Download schema"
        mkdir -p schema_local
        gsutil -m cp -r "~{schema_dir}*" schema_local/

        PT_ARG="";    [ -s "~{protein_table}" ] && PT_ARG="-t ~{protein_table}"
        GL_ARG="";    [ -s "~{genes_list}"    ] && GL_ARG="--gl ~{genes_list}"
        TAXA_ARG="";  [ -n "~{taxa}"          ] && TAXA_ARG="--taxa \"~{taxa}\""
        NS_ARG="";    [ "~{no_sparql}" = "true" ] && NS_ARG="--no-sparql"

        echo "[run] UniprotFinder"
        eval chewBBACA.py UniprotFinder \
            -g schema_local \
            -o annotations_output \
            --bsr ~{bsr} \
            --pm ~{proteome_matches} \
            --cpu ~{cpu} \
            $PT_ARG \
            $GL_ARG \
            $TAXA_ARG \
            $NS_ARG

        # Sortie principale
        cp annotations_output/schema_annotations.tsv schema_annotations.tsv || true
        tar -czf annotations_output.tar.gz annotations_output
    >>>

    output {
        File annotations_tsv = "schema_annotations.tsv"
        File annotations_tar = "annotations_output.tar.gz"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}

task compute_msa {
    input {
        String  schema_dir
        File    results_alleles
        File?   genes_list
        Int     translation_table = 11
        Int     cpu               = 8

        Int     disk_gb   = 50
        Int     memory_gb = 16
    }

    command <<<
        set -euo pipefail

        echo "[gsutil] Download schema"
        mkdir -p schema_local
        gsutil -m cp -r "~{schema_dir}*" schema_local/

        GL_ARG=""; [ -s "~{genes_list}" ] && GL_ARG="--gl ~{genes_list}"

        echo "[run] ComputeMSA"
        chewBBACA.py ComputeMSA \
            -g schema_local \
            -i "~{results_alleles}" \
            -o msa_output \
            --t ~{translation_table} \
            --cpu ~{cpu} \
            $GL_ARG

        tar -czf msa_output.tar.gz msa_output
    >>>

    output {
        File msa_tar = "msa_output.tar.gz"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}

task download_schema {
    input {
        String  species_id        
        String  schema_id        
        String  ns_url = "https://chewbbaca.online/api/NS/api/"
        File?   output_dir_name
        Int     disk_gb   = 50
        Int     memory_gb = 8
    }

    command <<<
        set -euo pipefail

        echo "[run] DownloadSchema (species=~{species_id}, schema=~{schema_id})"
        chewBBACA.py DownloadSchema \
            -sp "~{species_id}" \
            -sc "~{schema_id}" \
            -o downloaded_schema \
            --ns "~{ns_url}"

        tar -czf downloaded_schema.tar.gz downloaded_schema
    >>>

    output {
        File downloaded_schema_tar = "downloaded_schema.tar.gz"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        2
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 2
    }
}

task sync_schema {
    input {
        String  schema_dir          # GCP path vers schema local
        String  ns_url = "https://chewbbaca.online/api/NS/api/"
        String? submit              # "yes" to submit new alleles

        Int     disk_gb   = 50
        Int     memory_gb = 8
    }

    command <<<
        set -euo pipefail

        echo "[gsutil] Download schema"
        mkdir -p schema_local
        gsutil -m cp -r "~{schema_dir}*" schema_local/

        SUBMIT_ARG=""; [ -n "~{submit}" ] && SUBMIT_ARG="--submit ~{submit}"

        echo "[run] SyncSchema"
        chewBBACA.py SyncSchema \
            -g schema_local \
            --ns "~{ns_url}" \
            $SUBMIT_ARG

        echo "[compress] updated schema"
        tar -czf schema_synced.tar.gz schema_local
    >>>

    output {
        File schema_synced_tar = "schema_synced.tar.gz"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        2
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}

task ns_stats {
    input {
        String  ns_url = "https://chewbbaca.online/api/NS/api/"
        String? species_id

        Int     disk_gb   = 10
        Int     memory_gb = 4
    }

    command <<<
        set -euo pipefail

        SP_ARG=""; [ -n "~{species_id}" ] && SP_ARG="-sp ~{species_id}"

        echo "[run] NSStats"
        chewBBACA.py NSStats \
            --ns "~{ns_url}" \
            $SP_ARG \
            -o ns_stats_output

        tar -czf ns_stats_output.tar.gz ns_stats_output || true
    >>>

    output {
        File ns_stats_tar = "ns_stats_output.tar.gz"
    }

    runtime {
        docker:     "ummidock/chewbbaca:v3.3.10"
        cpu:        1
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 2
    }
}