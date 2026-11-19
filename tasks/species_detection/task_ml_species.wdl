version 1.0

# Species detection with the ML model (canonical k-mers + RandomForest) for ONE sample.
# The scripts (kmer_utils.py, predict_one.py) AND the model (strain_model.pkl) are inside the Docker
# image (/opt/species_prediction), see the Dockerfile. `ml_model` is only needed to try another model.
# The outputs have the same shape as the outputs of task_kraken2 so that they are interchangeable.

task ml_species_predict {
  input {
    File   read1
    File   read2
    File?  ml_model                       # optional: replaces the model built into the image
    String sample_id
    Float? threshold                      # empty = threshold stored in the model
    String docker = "bioinfomoh/species-prediction:0.1"
    Int    disk_gb = ceil(size(read1, "GB") + size(read2, "GB")) + 10
  }

  command <<<
    set -euo pipefail
    python3 /opt/species_prediction/predict_one.py \
      --r1 ~{read1} --r2 ~{read2} \
      --model ~{default="/opt/species_prediction/strain_model.pkl" ml_model} \
      --sample_id ~{sample_id} \
      ~{"--threshold " + threshold}
  >>>

  output {
    String version          = read_string("version.txt")
    File   report           = "report.tsv"
    String species_detected = read_string("species.txt")     # species, or UNKNOWN if confidence < threshold
    String top_taxa         = read_string("top_taxa.txt")    # 3 best species with their probability
    String genus_detected   = read_string("genus.txt")       # genus, "" if UNKNOWN (for task_species_match)
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}
