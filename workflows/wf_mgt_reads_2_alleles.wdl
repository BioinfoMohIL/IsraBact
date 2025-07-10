version 1.0

workflow mgt_reads_to_alleles {
    input {
        File    read1
        File    read2
        String samplename
        String  species
        String? serotype
        
        Int cpu             = 12   
        Int memory          = 30
        Int disk_size       = 250
        Float blastident    = 82
        Float hspident      = 0.95 
        Float locusnlimit   = 0.7
        String docker_image = "bioinfomoh/mgt_reads2alleles:1"
    }

    call reads_to_alleles {
        input: 
            read1 = read1,
            read2 = read2,
            samplename = samplename,
            species = species,
            serotype = serotype,
            cpu = cpu,
            memory = memory,
            disk_size = disk_size,
            docker_image = docker_image
        
    }

    output {
        Array[File] out = reads_to_alleles.out_file
    }
}

task reads_to_alleles {
    input {
        File read1
        File read2
        String samplename
        String species
        String? serotype

        Int cpu
        Int memory
        Int disk_size
        String docker_image
        
    }

    command <<<
        # Normalize user input to match the standard "Genus species" format:
        # first character uppercase, all remaining characters lowercase
        modif_species=$(echo "~{species}" | tr '[:upper:]' '[:lower:]')
        normalized_species="$(echo "${modif_species^}")"

        # If no ref alleles file provided:
        # The ref are in the /app/species_specific_files folder
        ref_dir_path="/app/MGT_reads2alleles/species_specific_files"
        
        case "$normalized_species" in
            "Bordetella pertussis")
                ref_path="${ref_dir_path}/pertussis/Pertussis_ref_alleles.fasta"
                ;;
            "Vibrio cholerae")
                ref_path="${ref_dir_path}/vibrio/Vibrio_ref_alleles.fasta"
                ;;
            "Salmonella enterica")
                ref_path="${ref_dir_path}/salmonella/salmonella_ref_alleles.fasta"
                ;;
            *)
                echo "❌ Unknown species: $normalized_species"
                exit 1
                ;;
        esac
                
        reads_to_alleles \
            -t ~{cpu} \
            -i ~{read1},~{read2} \
            --refalleles ${ref_path} \
            --species "${normalized_species}" \
            --strainid ~{samplename} \
            --kraken_db /app/db/minikraken_20171013_4GB \
            -o ./ \
            --genome_min 3500000 \
            --n50_min 15000 \
            --blastident 85 \
            --hspident 0.95 \
            --locusnlimit 0.7 \
            ~{if defined(serotype) then "--serotype " + serotype else "--no_serotyping"}

        # if [ -d ~{samplename} ]; then
        #     tar -czf out.tar.gz out
        # else
        #     echo "[ERROR] Output directory not found"
        #     exit 1
        # fi

    >>>

    output {
        Array[File] out_file = glob("/app/~{samplename}/*")

    }

    runtime {
        docker: docker_image
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk ~{disk_size} SSD"
    }


}