version 1.0

import "../../../tasks/utilities/file_handling/task_get_file_size.wdl" as task_get_file_size


workflow concatenate_illumina_lanes {
  input {
    String samplename
    File read1
    File? read2
  
  }
  
  call task_get_file_size.get_file_size {
    input:
      samplename = samplename, 
      read1 = read1,
      read2 = read2
  }

 

  output {
    Float read1_file_size_mb = get_file_size.fwd_file_size
    Float read2_file_size_mb = get_file_size.rev_file_size
    
  }
}