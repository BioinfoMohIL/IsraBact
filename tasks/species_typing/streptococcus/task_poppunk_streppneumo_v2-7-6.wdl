version 1.0


task poppunk {
  meta {
    description: "Using poppunk with GPS (Global Pneumococcal Sequencing project) database for Streptococcus pneumoniae typing"
  }

  input {
    File assembly
    String samplename
    
    File? GPS_db
    File? GPS_external_cluster
    String docker = "staphb/poppunk:2.7.6"
    Int disk_size = 100
    Int cpu = 12
    Int memory = 16
  }

  command <<<
    # get version information
    poppunk --version | sed 's/poppunk //' | tee VERSION
    
    GPS_DB_NAME=$(basename ~{GPS_db} | sed 's|.tar.gz||')
    
    echo "${GPS_DB_NAME}" > GPS_DB_NAME

    tar -xzvf ~{GPS_db}

    # create input TSV
    echo -e "~{samplename}\t~{assembly}" > ~{samplename}_poppunk_input.tsv

    mkdir -p ~{samplename}_poppunk
    poppunk_assign \
      --db ${GPS_DB_NAME} \
      --external-clustering ~{GPS_external_cluster} \
      --query ~{samplename}_poppunk_input.tsv \
      --output ~{samplename}_poppunk

    # parse output CSV for GPSC (Global Pneumococcal Sequence Cluster)
    if [ -f ~{samplename}_poppunk/~{samplename}_poppunk_external_clusters.csv ]; then
      cut -d ',' -f 2 ~{samplename}_poppunk/~{samplename}_poppunk_external_clusters.csv | tail -n 1 > GPSC.txt

      # if GPSC is "NA", overwrite with helpful message
      if [[ "$(cat GPSC.txt)" == "NA" ]]; then
        echo "Potential novel GPS Cluster identified, please email globalpneumoseq@gmail.com to have novel clusters added to the database and a GPSC cluster name assigned after you have checked for low level contamination which may contribute to biased accessory distances." >GPSC.txt
      fi
    else
      echo "poppunk failed" > GPSC.txt
    fi

  >>>

  output {
    String  poppunk_gps_cluster = read_string("GPSC.txt")
    File?   poppunk_gps_external_cluster_csv = "~{samplename}_poppunk/~{samplename}_poppunk_external_clusters.csv"
    File?   poppunk_dists_npy = "~{samplename}_poppunk/~{samplename}_poppunk.dists.npy"  
    File?   poppunk_dists_pkl = "~{samplename}_poppunk/~{samplename}_poppunk.dists.pkl"  
    File?   poppunk_h5 = "~{samplename}_poppunk/~{samplename}_poppunk.h5"  
    String  poppunk_version = read_string("VERSION")
    String  poppunk_docker = docker
    String  poppunk_GPS_db_version = read_string("GPS_DB_NAME")
  }
  
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
