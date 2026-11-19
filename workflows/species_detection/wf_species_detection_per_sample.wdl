version 1.0

import "../../tasks/metagenomics/task_kraken2.wdl" as task_kraken
import "../../tasks/species_detection/task_ml_species.wdl" as task_ml
import "../../tasks/utilities/task_species_match.wdl" as task_match

workflow wf_species_detection_per_sample {

  input {
    File   read1
    File   read2
    String sample_id

    Boolean use_ml = false        # false = kraken2 (défaut), true = modèle ML
    File?   ml_model              # optionnel: remplace le modèle intégré à l'image (strain_model.pkl)
    Float?  ml_threshold          # vide = seuil enregistré dans le modèle
    String  ml_docker = "bioinfomoh/species-prediction:0.1"
  }

  if (!use_ml) {
    call task_kraken.kraken2 {
      input:
        read1 = read1,
        read2 = read2,
        sample_id = sample_id
    }
  }

  if (use_ml) {
    call task_ml.ml_species_predict {
      input:
        read1 = read1,
        read2 = read2,
        ml_model = ml_model,          # vide = modèle intégré à l'image
        sample_id = sample_id,
        threshold = ml_threshold,
        docker = ml_docker
    }
  }

  # Même logique de match (préfixe du sample_id -> espèce attendue) pour kraken2 et pour le ML
  call task_match.species_match {
    input:
      sample_id = sample_id,
      detected_genus = select_first([kraken2.genus_detected, ml_species_predict.genus_detected]),
      detected_species = select_first([kraken2.species_detected, ml_species_predict.species_detected])
  }

  output {
    String spec_detect_ps_version          = select_first([kraken2.version, ml_species_predict.version])
    File   spec_detect_ps_report           = select_first([kraken2.report, ml_species_predict.report])
    String spec_detect_ps_species_detected = select_first([kraken2.species_detected, ml_species_predict.species_detected])
    String spec_detect_ps_top_taxa         = select_first([kraken2.top_taxa, ml_species_predict.top_taxa])
    String spec_detect_ps_expected_species = species_match.expected_species
    String spec_detect_ps_match            = species_match.match
    String spec_detect_ps_scanner          = if use_ml then "ML" else "kraken2"
  }
}
