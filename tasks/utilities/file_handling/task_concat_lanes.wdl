version 1.0


task cat_lanes {
  input {
    String samplename
    
    File read1_lane1
    File read1_lane2
    File? read1_lane3
    File? read1_lane4

    File? read2_lane1
    File? read2_lane2
    File? read2_lane3
    File? read2_lane4

    Int cpu = 2
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.2"
    Int memory = 4
  }
  meta {
    volatile: true
  }
  command <<<
    # exit task if anything throws an error (important for proper gzip format)
    set -euo pipefail
    
    exists() { [[ -f $1 ]]; }

    cat ~{read1_lane1} ~{read1_lane2} ~{read1_lane3} ~{read1_lane4} > "~{samplename}_merged_R1.fastq.gz"

    # Fetch the file size in MB
    fwd_file_size=$(stat -c%s "~{samplename}_merged_R1.fastq.gz")
    fwd_file_size_mb=$(awk -v size="$fwd_file_size" 'BEGIN {printf "%.2f", size / (1024*1024)}')
    echo "$fwd_file_size_mb" > fwd_size.txt
        
    if exists "~{read2_lane1}" ; then
      cat ~{read2_lane1} ~{read2_lane2} ~{read2_lane3} ~{read2_lane4} > "~{samplename}_merged_R2.fastq.gz"
    
      # Fetch the file size in MB
      rev_file_size=$(stat -c%s "~{samplename}_merged_R2.fastq.gz")
      rev_file_size_mb=$(awk -v size="$rev_file_size" 'BEGIN {printf "%.2f", size / (1024*1024)}')
      echo "$rev_file_size_mb" > rev_size.txt
    fi

    # ensure newly merged FASTQs are valid gzipped format
    gzip -t *merged*.gz
  >>>
  
  output {
    File read1_concatenated = "~{samplename}_R1.fastq.gz"
    File? read2_concatenated = "~{samplename}_R2.fastq.gz"
    
    Float fwd_file_size = read_float("fwd_size.txt")
    Float rev_file_size = read_float("rev_size.txt")
  }

  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}

