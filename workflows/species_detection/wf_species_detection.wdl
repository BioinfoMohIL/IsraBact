version 1.0

import "wf_species_detection_per_sample.wdl" as per_sample
import "../../tasks/utilities/task_plate_check.wdl" as task_plate

# Species detection for a SET of samples (run on a sample_set). It runs wf_species_detection_per_sample
# on every sample and, if a sending file is given, builds the plate check xlsx once all detections are done.
#
# In Terra (root entity = the set), for a table called "seq_data" for example:
#   read1       = this.seq_datas.read1
#   read2       = this.seq_datas.read2
#   sample_ids  = this.seq_datas.seq_data_id
#   sending_file = "gs://.../sending_file.xlsx"     (optional)
#
# To keep one column per sample in the data table, keep running wf_species_detection_per_sample on the
# samples themselves: this workflow returns its results as arrays on the set.

workflow wf_species_detection {

  input {
    Array[File]   read1
    Array[File]   read2
    Array[String] sample_ids

    File?   sending_file          # if set, the plate check xlsx is produced after the detection

    Boolean use_ml = false
    File?   ml_model
    Float?  ml_threshold
    String  ml_docker = "bioinfomoh/species-prediction:1.0"
  }

  scatter (i in range(length(sample_ids))) {
    call per_sample.wf_species_detection_per_sample as detect {
      input:
        read1 = read1[i],
        read2 = read2[i],
        sample_id = sample_ids[i],
        use_ml = use_ml,
        ml_model = ml_model,
        ml_threshold = ml_threshold,
        ml_docker = ml_docker
    }
  }

  if (defined(sending_file)) {
    call task_plate.plate_check {
      input:
        sending_file = select_first([sending_file]),
        sample_ids = sample_ids,
        species_detected = detect.sp_detection_species_detected
    }
  }

  output {
    Array[String] sp_detection_scanner          = detect.sp_detection_scanner
    Array[File]   sp_detection_report           = detect.sp_detection_report
    File?         sp_detection_plate_check      = plate_check.report_xlsx
  }
}
