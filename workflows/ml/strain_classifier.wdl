version 1.0

# Entraînement d'un classifieur d'espèce (k-mers canoniques + RandomForest).
#   1. scatter  : une tâche `featurize` par échantillon (lecture des FASTQ, petite VM)
#   2. gather   : une tâche `train` (CV out-of-fold + modèle final) sur les petits vecteurs
#
# Dans Terra: lancer sur un sample_set, avec par exemple
#   read1   = this.samples.read1
#   read2   = this.samples.read2
#   species = this.samples.species   (optionnel: sans ça, l'espèce est déduite du préfixe)

workflow strain_classifier_train {
  input {
    Array[File]   read1
    Array[File]   read2
    Array[String] species = []   # optionnel: si vide, déduit du préfixe du nom d'échantillon

    File kmer_utils_py
    File featurize_py
    File train_py

    Int    k           = 8
    Int    max_reads   = 100000
    Int    folds       = 5
    Float  threshold   = 0.5
    File?  groups_csv            # optionnel: CSV sample,group (run/batch)
    String docker = "bioinfomoh/strain-clf:0.1"
  }

  scatter (i in range(length(read1))) {
    String sample_id = sub(basename(read1[i]), "(_R1.*|_1\\.f.*)$", "")

    call featurize {
      input:
        read1         = read1[i],
        read2         = read2[i],
        kmer_utils_py = kmer_utils_py,
        featurize_py  = featurize_py,
        k             = k,
        max_reads     = max_reads,
        docker        = docker
    }
  }

  call train {
    input:
      sample_ids    = sample_id,
      species       = species,
      features      = featurize.features,
      n_reads       = featurize.n_reads,
      kmer_utils_py = kmer_utils_py,
      train_py      = train_py,
      k             = k,
      max_reads     = max_reads,
      folds         = folds,
      threshold     = threshold,
      groups_csv    = groups_csv,
      docker        = docker
  }

  output {
    File model           = train.model
    File oof_predictions = train.oof_predictions
    File cv_report       = train.cv_report
  }
}

task featurize {
  input {
    File read1
    File read2
    File kmer_utils_py
    File featurize_py
    Int  k
    Int  max_reads
    String docker
    Int  disk_gb = ceil(size(read1, "GB") + size(read2, "GB")) + 10
  }

  command <<<
    set -euo pipefail
    cp ~{kmer_utils_py} kmer_utils.py
    cp ~{featurize_py} featurize.py
    python3 featurize.py \
      --r1 ~{read1} --r2 ~{read2} \
      --k ~{k} --max_reads ~{max_reads} \
      --out features.npy --n_reads_out n_reads.txt
  >>>

  output {
    File features = "features.npy"
    Int  n_reads  = read_int("n_reads.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task train {
  input {
    Array[String] sample_ids
    Array[String] species = []
    Array[File]   features
    Array[Int]    n_reads
    File   kmer_utils_py
    File   train_py
    Int    k
    Int    max_reads
    Int    folds
    Float  threshold
    File?  groups_csv
    String docker
  }

  command <<<
    set -euo pipefail
    cp ~{kmer_utils_py} kmer_utils.py
    cp ~{train_py} train.py

    # manifest: une ligne par échantillon (les chemins de `features` sont ceux localisés)
    { printf 'sample\tspecies\tfeatures_path\tn_reads\n'; \
      paste ~{write_lines(sample_ids)} ~{write_lines(species)} \
            ~{write_lines(features)} ~{write_lines(n_reads)}; } > manifest.tsv

    python3 train.py \
      --features_tsv manifest.tsv --output out \
      --k ~{k} --max_reads ~{max_reads} --folds ~{folds} --threshold ~{threshold} \
      ~{"--groups_csv " + groups_csv} \
      | tee cv_report.txt
  >>>

  output {
    File model           = "out/strain_model.pkl"
    File oof_predictions = "out/oof_predictions.csv"
    File cv_report       = "cv_report.txt"
  }

  runtime {
    docker: docker
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 20 HDD"
    preemptible: 1
  }
}
