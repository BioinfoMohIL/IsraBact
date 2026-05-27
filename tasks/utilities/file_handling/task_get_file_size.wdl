version 1.0


task get_file_size {
  input {    
    String samplename
    File read1
    File? read2
    Int disk_size = 50
   
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.2"
  }
  meta {
    volatile: true
  }
  command <<<
    # exit task if anything throws an error (important for proper gzip format)
    set -euo pipefail
    
    exists() { [[ -f $1 ]]; }

    cat ~{read1} > "~{samplename}_R1.fastq.gz"

    # Fetch the file size in MB
    fwd_file_size=$(stat -c%s "~{samplename}_R1.fastq.gz")
    fwd_file_size_mb=$(awk -v size="$fwd_file_size" 'BEGIN {printf "%.2f", size / (1024*1024)}')
    echo "$fwd_file_size_mb" > fwd_size.txt
        
    if exists "~{read2}" ; then
      cat ~{read2} > "~{samplename}_R2.fastq.gz"
    
      # Fetch the file size in MB
      rev_file_size=$(stat -c%s "~{samplename}_R2.fastq.gz")
      rev_file_size_mb=$(awk -v size="$rev_file_size" 'BEGIN {printf "%.2f", size / (1024*1024)}')
      echo "$rev_file_size_mb" > rev_size.txt
    fi

  >>>
  
  output {
    Float fwd_file_size = read_float("fwd_size.txt")
    Float rev_file_size = read_float("rev_size.txt")
  }

  runtime {
    docker: "~{docker}"
    disks: "local-disk " + disk_size + " SSD"
    preemptible: 1
  }
}

