version 1.0

task shovill {
  input {
    File read1
    File read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/shovill:1.1.0"
    Int disk_size = 100
    Int cpu = 16
    Int memory = 16
    Int? depth
    String? genome_length
    Int min_contig_length = 200
    Float? min_coverage
    String assembler = "skesa"
    String? assembler_options
    String? kmers
    Boolean trim = false
    Boolean noreadcorr = false
    Boolean nostitch = false
    Boolean nocorr = false
  }

  command <<<
    shovill --version | head -1 | tee VERSION
    shovill \
      --outdir out \
      --R1 ~{read1} \
      --R2 ~{read2} \
      --minlen ~{min_contig_length} \
      ~{'--depth ' + depth} \
      ~{'--gsize ' + genome_length} \
      ~{'--mincov ' + min_coverage} \
      ~{'--assembler ' + assembler} \
      ~{'--opts ' + assembler_options} \
      ~{'--kmers ' + kmers} \
      ~{true='--trim' false='' trim} \
      ~{true='--noreadcorr' false='' noreadcorr} \
      ~{true='--nostitch' false='' nostitch} \
      ~{true='--nocorr' false='' nocorr} \
      --cpus ~{cpu}

    mv out/contigs.fa out/~{samplename}_contigs.fasta

    if [ "~{assembler}" == "spades" ] ; then
      mv out/contigs.gfa out/~{samplename}_contigs.gfa
    elif [ "~{assembler}" == "megahit" ] ; then
      mv out/contigs.fastg out/~{samplename}_contigs.fastg
    elif [ "~{assembler}" == "velvet" ] ; then
      mv out/contigs.LastGraph out/~{samplename}_contigs.LastGraph
    fi
  >>>
  output {
    File assembly_fasta = "out/~{samplename}_contigs.fasta"
    File? contigs_gfa = "out/~{samplename}_contigs.gfa"
    File? contigs_fastg = "out/~{samplename}_contigs.fastg"
    File? contigs_lastgraph = "out/~{samplename}_contigs.LastGraph"
    String shovill_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
