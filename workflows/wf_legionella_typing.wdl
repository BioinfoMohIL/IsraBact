version 1.0

workflow legionella_typing {
    meta {
        description: "SBT Legionella pneumophila isolates."
        author: "David Maimoun"
        organization: "MOH Jerusalem"
    }

    input {
        String samplename
        File read1 
        File read2
        String docker = "staphb/elgato:1.21.2"
    }

    call elgato_reads {
        input:
        read1       = read1,
        read2       = read2,
        samplename  = samplename,
        docker      = docker
    }

    output {
        String sbt_elgato_version = elgato_reads.elgato_version
        String sbt_elgato_reads = elgato_reads.sbt
        File sbt_elgato_possible_sts = elgato_reads.possible_mlsts
        File sbt_elgato_inter_out = elgato_reads.intermediate_outputs
        File sbt_elgato_alleles = elgato_reads.alleles
    }
}

task elgato_reads {
  input {
    File read1
    File read2
    String samplename
    String docker
    
    Int cpu = 8
  }

  command <<<
    el_gato.py -v > VERSION
    
    el_gato.py --read1 ~{read1} --read2 ~{read2} --out ./out
  
    st=$(awk -F "\t" 'NR==2 {print $2}' ./out/possible_mlsts.txt)
    
    if [ -z "$st" ]; then
      st="No ST predicted!"
    else
      st="ST"$st
    fi
    
    echo $st > SBT

    mv out/possible_mlsts.txt ~{samplename}_possible_mlsts.txt
    mv out/intermediate_outputs.txt ~{samplename}_intermediate_outputs.txt
    mv out/identified_alleles.fna ~{samplename}_identified_alleles.fna
  >>>

  output {    
    String elgato_version = read_string("VERSION")
    String sbt = read_string("SBT")
    File possible_mlsts = "~{samplename}_possible_mlsts.txt"
    File intermediate_outputs = "~{samplename}_intermediate_outputs.txt"
    File alleles = "~{samplename}_identified_alleles.fna"
  }

  runtime {
    docker: docker
    memory: "~{cpu} GB"
    cpu: cpu
  }

}





