version 1.0

# Compares the detected species of all samples with the sending file (plate / well / expected species)
# and writes an xlsx with one sheet per plate that has results: an 8x12 grid of sample names coloured
#   green = match, orange = same genus (different species), red = mismatch, grey = no result,
# with the list of errors below the grid. Handy to spot sample mix-ups made during sequencing.
#
# plate_check.py (standard library only) is built into the Docker image: /opt/species_prediction/plate_check.py

task plate_check {
  input {
    File          sending_file         # the sending file (.xlsx); the table can start on any row
    Array[String] sample_ids           # same order as species_detected
    Array[String] species_detected
    String        output_name = "plate_check.xlsx"
    String        docker = "bioinfomoh/species-prediction:0.1"
  }

  command <<<
    set -euo pipefail
    { printf 'sample_id\tspecies_detected\n'; \
      paste ~{write_lines(sample_ids)} ~{write_lines(species_detected)}; } > results.tsv

    python3 /opt/species_prediction/plate_check.py \
      --sending "~{sending_file}" \
      --results results.tsv --id-col sample_id --species-col species_detected \
      --out "~{output_name}"
  >>>

  output {
    File report_xlsx = "~{output_name}"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}
