version 1.0

# Détection d'espèce par le modèle ML (k-mers canoniques + RandomForest) pour UN échantillon.
# Les scripts (kmer_utils.py, predict_one.py) ET le modèle (strain_model.pkl) sont dans l'image Docker
# (/opt/species_prediction), voir Dockerfile. `ml_model` n'est utile que pour tester un autre modèle.
# Les sorties ont la même forme que celles de task_kraken2 pour pouvoir être interchangées.

task ml_species_predict {
  input {
    File   read1
    File   read2
    File?  ml_model                       # optionnel: remplace le modèle intégré à l'image
    String sample_id
    Float? threshold                      # vide = seuil enregistré dans le modèle
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
    String species_detected = read_string("species.txt")     # espèce, ou UNKNOWN si confiance < seuil
    String top_taxa         = read_string("top_taxa.txt")    # 3 meilleures espèces avec leur probabilité
    String genus_detected   = read_string("genus.txt")       # genre, "" si UNKNOWN (pour task_species_match)
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}
