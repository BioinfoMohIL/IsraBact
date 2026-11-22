version 1.0

import "wf_species_detection_per_sample.wdl" as per_sample
import "../../tasks/utilities/task_species_report.wdl" as task_report
import "../../tasks/utilities/task_plate_check.wdl" as task_plate

# Species detection for a SET of samples (run on a sample_set). It runs wf_species_detection_per_sample
# on every sample, then returns ONE report for the whole run plus a one-line summary, instead of
# long per-sample arrays. If a sending file is given, it also builds the plate check xlsx.
#
# In Terra (root entity = the set), for a table called "seq_data" for example:
#   read1        = this.seq_datas.read1
#   read2        = this.seq_datas.read2
#   sample_ids   = this.seq_datas.seq_data_id
#   sending_file = "gs://.../sending_file.xlsx"     (optional)
#
# Outputs (set level):
#   sp_detection_scanner      "kraken2" or "ML"
#   sp_detection_report       one CSV: sample_id,detected_species,match,expected_species,scanner,top_taxa
#   sp_detection_summary      "ALL MATCHED!" or e.g. "MATCH: 90,APPROX: 2,NOMATCH: 3,UNKNOWN: 1"
#   sp_detection_plate_check  xlsx, one sheet per plate (only if sending_file is set; empty if it failed)
#   sp_detection_plate_status "OK" or "FAILED: reason" (empty if no sending_file). A failed plate check
#                             does not fail the workflow: use the report with plate_check_app.py instead.
#
# To keep one column per sample in the data table, run wf_species_detection_per_sample on the samples.

workflow wf_species_detection {

  input {
    Array[File]   read1
    Array[File]   read2
    Array[String] sample_ids

    File?   sending_file          # if set, the plate check xlsx is produced after the detection

    Boolean use_ml = false
    File?   ml_model
    Float?  ml_threshold
    String  ml_docker = "bioinfomoh/species-detection:2.0"
  }

  String scanner = if use_ml then "ML" else "kraken2"

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

  call task_report.species_report {
    input:
      sample_ids = sample_ids,
      species_detected = detect.sp_detection_species_detected,
      matches = detect.sp_detection_match,
      expected_species = detect.sp_detection_expected_species,
      top_taxa = detect.sp_detection_top_taxa,
      scanner = scanner
  }

  call task_report.summarize_species {
    input:
      report_csv = species_report.report
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
    String sp_detection_scanner     = scanner
    File   sp_detection_report      = species_report.report
    String sp_detection_summary     = summarize_species.summary
    File?  sp_detection_plate_check = plate_check.report_xlsx
    String? sp_detection_plate_status = plate_check.status
  }
}
