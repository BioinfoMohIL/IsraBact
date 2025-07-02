version 1.0

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

    mkdir out
    fq2dna | grep "fq2dna v" | awk '{print $3}' | tr -d 'v' | tee "VERSION"
    
    fq2dna \
        -1 ~{read1} \
        -2 ~{read2} \
        -o out \
        -b ${basename} \
        -s ~{organism} \
        -T "${species_formatted}" \
        -a ~{alien_tag}

    tar -czvf txt_info_files.tar.gz out/*.txt

    
  >>>

  output {
    String version           = read_string('VERSION')
    File assembly_fasta      = "out/-s.all.fasta"
    File metrics_zip         = "txt_info_files.tar.gz"
    File selected_scaffolds  = "out/-s.scf.fasta"
    File selected_contigs    = "out/-s.agp.fasta"
    File scaffolding_info    = "out/-s.agp"
  }

  runtime {
    docker: "bioinfomoh/fq2dna:25.03"
    cpu: 20
    memory: "16G"
    disks: "local-disk 10 HDD"
  }
}
