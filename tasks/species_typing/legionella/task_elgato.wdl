version 1.0

task elgato_reads {
  input {
    File read1
    File read2
    String samplename
    String docker = "staphb/elgato:1.21.2"
    
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
    String elgato_sbt = read_string("SBT")
    File elgato_possible_mlsts = "~{samplename}_possible_mlsts.txt"
    File elgato_intermediate_outputs = "~{samplename}_intermediate_outputs.txt"
    File elgato_alleles = "~{samplename}_identified_alleles.fna"
  }

  runtime {
    docker: docker
    memory: "~{cpu} GB"
    cpu: cpu
  }

}





