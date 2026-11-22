version 1.0

# Compares the detected species of all samples with the sending file (plate / well / expected species)
# and writes an xlsx with one sheet per plate that has results: an 8x12 grid of sample names coloured
#   green = match, orange = same genus (different species), red = mismatch, grey = no result,
# with the list of errors below the grid. Handy to spot sample mix-ups made during sequencing.
#
# plate_check.py (standard library only) is built into the Docker image: /opt/species_prediction/plate_check.py
#
# This step is a "nice to have": if plate_check.py fails (bad sending file, unexpected input...), the task
# does NOT fail, so the rest of the workflow (report, summary) is kept. It only prints a warning, and
# `status` says what happened ("OK" or "FAILED: ..."). The xlsx is then missing (File? = not set):
# use the Streamlit app (plate_check_app.py) with the report instead.
# (A failure of the VM itself, e.g. the Docker image cannot be pulled, cannot be caught here.)

task plate_check {
  input {
    File          sending_file         # the sending file (.xlsx); the table can start on any row
    Array[String] sample_ids           # same order as species_detected
    Array[String] species_detected
    String        output_name = "plate_check.xlsx"
    File?         plate_check_script   # optional: plate_check.py from a bucket (any image with python3 then works)
    String        docker = "ghcr.io/bioinfomoh/species_detection:2.0"
  }

  command <<<
    set -uo pipefail
    { printf 'sample_id\tspecies_detected\n'; \
      paste ~{write_lines(sample_ids)} ~{write_lines(species_detected)}; } > results.tsv

    if python3 ~{default="/opt/species_prediction/plate_check.py" plate_check_script} \
         --sending "~{sending_file}" \
         --results results.tsv --id-col sample_id --species-col species_detected \
         --out "~{output_name}" 2> err.txt; then
        echo "OK" > status.txt
    else
        cat err.txt >&2
        echo "WARNING: plate check failed, the workflow continues without it" >&2
        rm -f "~{output_name}"
        echo "FAILED: $(tail -n 1 err.txt | tr -d '\n' | cut -c1-200)" > status.txt
    fi
    exit 0
  >>>

  output {
    File?  report_xlsx = "~{output_name}"     # not set if the plate check failed
    String status      = read_string("status.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 2
    maxRetries: 1            # retries a VM / image-pull failure once
  }
}
