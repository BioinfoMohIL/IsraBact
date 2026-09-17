task cat_files {
  input {
    Array[File] files_to_cat
    String concatenated_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Boolean get_samplename
    Boolean skip_extra_headers = false
    String delimiter = "\t"
  }

  meta {
    # added so that call caching is always turned off
    volatile: true
  }

  command <<<
    file_array=(~{sep=' ' files_to_cat})
    touch ~{concatenated_file_name}

    for index in ${!file_array[@]}; do
      file=${file_array[$index]}

      if ~{get_samplename}; then
        # extrait le samplename depuis le nom de fichier (ex: NM0001_results.tsv -> NM0001)
        samplename=$(basename "${file}" | cut -d'_' -f1)

        if [ "$index" -eq 0 ]; then
          if ~{skip_extra_headers}; then
            # 1er fichier : on garde le header, renomme en "samplename", et ajoute la colonne aux autres lignes
            awk -v var="$samplename" 'BEGIN{FS=OFS="~{delimiter}"} NR==1{print "samplename", $0; next} {print var, $0}' "${file}" >> ~{concatenated_file_name}
          else
            awk -v var="$samplename" 'BEGIN{FS=OFS="~{delimiter}"} {print var, $0}' "${file}" >> ~{concatenated_file_name}
          fi
        else
          if ~{skip_extra_headers}; then
            tail -n +2 "${file}" | awk -v var="$samplename" 'BEGIN{FS=OFS="~{delimiter}"} {print var, $0}' >> ~{concatenated_file_name}
          else
            awk -v var="$samplename" 'BEGIN{FS=OFS="~{delimiter}"} {print var, $0}' "${file}" >> ~{concatenated_file_name}
          fi
        fi
      else
        # comportement d'origine, sans colonne samplename
        if ! ~{skip_extra_headers} ; then
          cat ${file} >> ~{concatenated_file_name}
        else
          if [ $index == 0 ]; then
            cat ${file} >> ~{concatenated_file_name}
          else
            tail -n +2 ${file} >> ~{concatenated_file_name}
          fi
        fi
      fi
    done
  >>>
  output {
    File concatenated_files = "~{concatenated_file_name}"
  }
  runtime {
    docker: "~{docker_image}"
    memory:  "8 GB"
    cpu: 2
    disks: "local-disk 100 HDD"
    disk: "100 GB"
    preemptible: 0
  }
}