version 1.0

task fastqc {
  input {
    File read1
    File read2

    Int memory = 4
    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/fastqc:0.12.1"
  }
  String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
  String read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
  command <<<
      # get fastqc version
    fastqc --version | tee VERSION
    
    # run fastqc: 
    # --extract: uncompress output files
    fastqc \
      --outdir . \
      --threads ~{cpu} \
      --extract \
      ~{read1} \
      ~{read2}

    grep "Total Sequences" ~{read1_name}_fastqc/fastqc_data.txt | cut -f 2 | tee READ1_SEQS
    read1_seqs=$(cat READ1_SEQS)
    grep "Total Sequences" ~{read2_name}_fastqc/fastqc_data.txt | cut -f 2 | tee READ2_SEQS
    read2_seqs=$(cat READ2_SEQS)

    # capture number of read pairs
    if [ "${read1_seqs}" == "${read2_seqs}" ]; then
      read_pairs=${read1_seqs}
    else
      read_pairs="Uneven pairs: R1=${read1_seqs}, R2=${read2_seqs}"
    fi
    
    echo "$read_pairs" | tee READ_PAIRS
  >>>

  output {
    File read1_fastqc_html = "~{read1_name}_fastqc.html"
    File read1_fastqc_zip = "~{read1_name}_fastqc.zip"
    File read2_fastqc_html = "~{read2_name}_fastqc.html"
    File read2_fastqc_zip = "~{read2_name}_fastqc.zip"
    
    Int read1_seq = read_int("READ1_SEQS")
    Int read2_seq = read_int("READ2_SEQS")
    String read_pairs = read_string("READ_PAIRS")
    String version = read_string("VERSION")
    String fastqc_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}

task busco {
  meta {
    description: "Run BUSCO on assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/ezlabgva/busco:v5.7.1_cv1"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
    Boolean eukaryote = false
  }
  command <<<
    # get version
    busco --version | tee "VERSION"
 
    # run busco
    # -i input assembly
    # -m geno for genome input
    # -o output file tag
    # --auto-lineage-euk looks at only eukaryotic organisms
    # --auto-lineage-prok looks at only prokaryotic organisms; default
    busco \
      -i ~{assembly} \
      -c ~{cpu} \
      -m geno \
      -o ~{samplename} \
      ~{true='--auto-lineage-euk' false='--auto-lineage-prok' eukaryote}

    # check for existence of output file; otherwise display a string that says the output was not created
    if [ -f ~{samplename}/short_summary.specific.*.~{samplename}.txt ]; then

      # grab the database version and format it according to BUSCO recommendations
      # pull line out of final specific summary file
      # cut out the database name and date it was created
      # sed is to remove extra comma and to add parentheses around the date and remove all tabs
      # finally write to a file called DATABASE
      cat ~{samplename}/short_summary.specific.*.~{samplename}.txt | grep "dataset is:" | cut -d' ' -f 6,9 | sed 's/,//; s/ / (/; s/$/)/; s|[\t]||g' | tee DATABASE
      
      # extract the results string; strip off all tab and space characters; write to a file called BUSCO_RESULTS
      cat ~{samplename}/short_summary.specific.*.~{samplename}.txt | grep "C:" | sed 's|[\t]||g; s| ||g' | tee BUSCO_RESULTS

      # rename final output file to predictable name
      cp -v ~{samplename}/short_summary.specific.*.~{samplename}.txt ~{samplename}_busco-summary.txt
    else
      echo "BUSCO FAILED" | tee BUSCO_RESULTS
      echo "NA" > DATABASE
    fi
  >>>
  output {
    String busco_version = read_string("VERSION")
    String busco_docker = docker
    String busco_database = read_string("DATABASE")
    String busco_results = read_string("BUSCO_RESULTS")
    File   busco_report = "~{samplename}_busco-summary.txt"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}
