version 1.0

task meningotype {
  meta {
    description: "Serotyping, finetyping, BAST and MenDeVAR analysis of Neisseria meningitidis"
    version: "1.0"
  }

  input {
    File assembly
    String samplename

    String docker = "bioinfomoh/meningotype:1.0"

    Int disk_size = 100
    Int cpu = 8
    Int memory = 8
  }

  command <<<
    set -euo pipefail

    # Tool version
    meningotype --version 2>&1 | tee VERSION

    # Run complete Neisseria typing
    meningotype \
      --all \
      --cpus ~{cpu} \
      ~{assembly} \
      > ~{samplename}.tsv

    # Extract results
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $2}'  > MENINGOTYPE_SEROTYPE
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $4}'  > MENINGOTYPE_MLST
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $5}'  > MENINGOTYPE_PORA
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $6}'  > MENINGOTYPE_FETA
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $7}'  > MENINGOTYPE_PORB
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $8}'  > MENINGOTYPE_FHBP
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $9}'  > MENINGOTYPE_NHBA
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $10}' > MENINGOTYPE_NADA
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $11}' > MENINGOTYPE_BAST

    # MenDeVAR results
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $12}' > MENINGOTYPE_BEXSERO_MENDEVAR
    tail -n 1 ~{samplename}.tsv | awk -F '\t' '{print $13}' > MENINGOTYPE_TRUMENBA_MENDEVAR

  >>>

  output {
    File meningotype_tsv = "~{samplename}.tsv"

    String meningotype_version = read_string("VERSION")

    String meningotype_serogroup = read_string("MENINGOTYPE_SEROTYPE")
    String meningotype_MLST = read_string("MENINGOTYPE_MLST")

    String meningotype_PorA = read_string("MENINGOTYPE_PORA")
    String meningotype_FetA = read_string("MENINGOTYPE_FETA")
    String meningotype_PorB = read_string("MENINGOTYPE_PORB")

    String meningotype_fHbp = read_string("MENINGOTYPE_FHBP")
    String meningotype_NHBA = read_string("MENINGOTYPE_NHBA")
    String meningotype_NadA = read_string("MENINGOTYPE_NADA")

    String meningotype_BAST = read_string("MENINGOTYPE_BAST")

    String meningotype_Bexsero_MenDeVAR = read_string("MENINGOTYPE_BEXSERO_MENDEVAR")
    String meningotype_Trumenba_MenDeVAR = read_string("MENINGOTYPE_TRUMENBA_MENDEVAR")
  }

  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
  }
}
