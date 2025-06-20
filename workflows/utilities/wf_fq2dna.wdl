version 1.0

workflow wf_fq2dna {
  input {
    File read1
    File read2
    String genus
    String species
    String strain 
    String? organism = "B" # "B(Bacteria), P(Prokariote),  E(Eukaryote), S(Standard), V(Virus)"
    String? alien_tag = "AUTO"
  }

  call fq2dna_run {
    input:
        read1 = read1,
        read2 = read2,
        genus = genus,
        species = species,
        strain = strain ,
        organism = organism,
        alien_tag = alien_tag
    }

#   output {
#     File fq2dna_report = fq2dna_run.report
#   }
}

task fq2dna_run {
  input {
    File read1
    File read2
    String genus
    String species
    String strain 
    String? organism
    String? alien_tag 

  }

  command <<< 
    genus="~{genus}"
    species="~{species}"
    strain="~{strain}"
    basename=$(echo "${genus:0:1}${species:0:1}${strain}" | tr '[:upper:]' '[:lower:]')
    species="${genus} ${species} ${strain}"
    
    fq2dna \
        -1 ~{read1} \
        -2 ~{read2} \
        -o out \
        -b ${basename} \
        -s ~{organism} \
        -T "${species}" \
        -a ~{alien_tag}

    echo out > report.txt
  >>>

  output {
    File report = "report.txt"
  }

  runtime {
    docker: "bioinfomoh/fq2dna:1"
    cpu: 8
    memory: "8G"
    disks: "local-disk 10 HDD"
  }
}
