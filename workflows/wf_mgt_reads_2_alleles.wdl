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

        Int min_largest_contig 
        Int max_contig_no 
        Int genome_min   
        Int n50_min         
        Int blastident   
        Float hspident          
        Float locusnlimit 
          
        # Int min_largest_contig = 50000
        # Int max_contig_no = 700
        # Int genome_min    = 3500000
        # Int n50_min       = 15000    
        # Int blastident    = 85
        # Float hspident      = 0.95     
        # Float locusnlimit   = 0.7
 
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
            docker_image = docker_image,

            min_largest_contig  = min_largest_contig,
            max_contig_no = max_contig_no,
            genome_min = genome_min,
            n50_min = n50_min,   
            blastident = blastident, 
            hspident = hspident,     
            locusnlimit = locusnlimit
        
    }

    output {
        File mgt_alleles = reads_to_alleles.alleles
        File mgt_contigs = reads_to_alleles.contigs
        File mgt_pass = reads_to_alleles.pass
        File mgt_asssembly_stats = reads_to_alleles.asssembly_stats
        File mgt_kraken_report = reads_to_alleles.kraken_report
        File? mgt_shovill = reads_to_alleles.shovill
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

        Int genome_min
        Int min_largest_contig
        Int max_contig_no
        Int n50_min      
        Int blastident   
        Float hspident     
        Float locusnlimit  
        
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
            ~{'--min_largest_contig '  + min_largest_contig} \
            ~{' --max_contig_no '  +  max_contig_no} \
            ~{'--genome_min '  + genome_min} \
            ~{'--genome_min '  + genome_min} \
            ~{'--n50_min '     + n50_min} \
            ~{'--blastident '  + blastident} \
            ~{'--hspident '    + hspident } \
            ~{'--locusnlimit ' + locusnlimit } \
            ~{if defined(serotype) then "--serotype " + serotype else "--no_serotyping"}


        found_path=$(find / -type d -name "~{samplename}")
        if [[ -n "$found_path" ]]; then
            echo "found in ${found_path}"
            mv "$found_path" .
        else
            echo "❌ Directory '~{samplename}' not found."
        fi

        if [[ -d "~{samplename}/~{samplename}_shovill" ]]; then
            tar -czf "~{samplename}_shovill.tar.gz" "~{samplename}/~{samplename}_shovill"
            echo "✅ Created: ~{samplename}_shovill.tar.gz"
        else
            echo "❌ Directory not found: ~{samplename}/~{samplename}_shovill"
        fi


    >>>

    
    output {
       File alleles = "~{samplename}/~{samplename}_alleles.fasta"
       File contigs = "~{samplename}/~{samplename}_contigs.fa"
       File pass = "~{samplename}/~{samplename}_pass.fasta"
       File asssembly_stats = "~{samplename}/~{samplename}_assembly_stats.txt"
       File kraken_report = "~{samplename}/~{samplename}_kraken_out_report.txt"
       File? shovill = "~{samplename}_shovill.tar.gz"
    }

    runtime {
        docker: docker_image
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk ~{disk_size} SSD"
    }


}