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

#   output {
#     File fq2dna_report = fq2dna_run.report
#   }
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
  >>>

#   output {
#     File report = "out/report.txt"
#   }

  runtime {
    docker: "bioinfomoh/fq2dna:1"
    cpu: 20
    memory: "16G"
    disks: "local-disk 10 HDD"
  }
}
