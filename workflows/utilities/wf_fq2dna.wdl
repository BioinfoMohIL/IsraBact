version 1.0

workflow wf_fq2dna {
  input {
    File read1
    File read2
    String species
    String? organism = "B" # "B(Bacteria), P(Prokariote),  E(Eukaryote), S(Standard), V(Virus)"
    String? alien_tag = "AUTO"
  }

  call fq2dna_run {
    input:
        read1 = read1,
        read2 = read2,
        species = species,
        organism = organism,
        alien_tag = alien_tag
    }

  output {
    File fq2dna_assembly_fasta      = fq2dna_run.assembly_fasta
    File fq2dna_metrics_zip         = fq2dna_run.metrics_zip
    File fq2dna_selected_scaffolds  = fq2dna_run.selected_scaffolds
    File fq2dna_selected_contigs    = fq2dna_run.selected_contigs
    File fq2dna_scaffolding_info:   = fq2dna_run.scaffolding_info
  }
}


task fq2dna_run {
  input {
    File read1
    File read2
    String species
    String? organism
    String? alien_tag 

  }

  command <<< 
    genus=$(echo "~{species}" | awk '{print $1}')
    species_name=$(echo "~{species}" | awk '{print $2}')
    strain=$(echo "~{species}" | awk '{print $3}')

    species_formatted=$(echo "${genus:0:1}${species_name:0:4}${strain:0:2}" | tr '[:upper:]' '[:lower:]')

    fq2dna \
        -1 ~{read1} \
        -2 ~{read2} \
        -o out \
        -b ${basename} \
        -s ~{organism} \
        -T "${species_formatted}" \
        -a ~{alien_tag}

    zip txt_info_files.zip out/*.txt
  >>>

  output {
    File assembly_fasta      = "out/-s.all.fasta"
    File metrics_zip         = "txt_info_files.zip"
    File selected_scaffolds  = "out/-s.scf.fasta"
    File selected_contigs    = "out/-s.agp.fasta"
    File scaffolding_info:   = "out/-s.agp"
  }

  runtime {
    docker: "bioinfomoh/fq2dna:1"
    cpu: 20
    memory: "16G"
    disks: "local-disk 10 HDD"
  }
}
