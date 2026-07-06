version 1.0

import "../../tasks/metagenomics/task_kraken2.wdl" as task_kraken
import "../../tasks/utilities/file_handling/decompress_file.wdl" as task_unzip 

workflow wf_kraken2 {
  
  input {
    File read1
    File read2
    String sample_id
    File kraken_db = "gs://fc-5d4556f8-3de6-4709-85da-11445772644d/db/minikraken2_v2_8GB_201904_UPDATE.zip"
  }

  call task_unzip.decompress {
    input:
        archive_file = kraken_db
  }

  call task_kraken.kraken2 {
    input:
        read1 = read1,
        read2 = read2,
        sample_id = sample_id,
        kraken_db_path = decompress.decompressed_dir_path
  }

  output {
        String version = kraken2.version
        File report = kraken2.report
        String species_detected = kraken2.species_detected
    }
}
  