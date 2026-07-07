version 1.0

import "../../tasks/metagenomics/task_kraken2.wdl" as task_kraken

workflow wf_kraken2 {
  
  input {
    File read1
    File read2
    String sample_id
  }

  call task_kraken.kraken2 {
    input:
        read1 = read1,
        read2 = read2,
        sample_id = sample_id
  }

  output {
        String version = kraken2.version
        File report = kraken2.report
        String species_detected = kraken2.species_detected
    }
}
  