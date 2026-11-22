version 1.0

import "../../tasks/metagenomics/task_kraken2.wdl" as task_kraken
import "../../tasks/species_detection/task_ml_species.wdl" as task_ml
import "../../tasks/utilities/task_species_match.wdl" as task_match

# Species detection for ONE sample, with kraken2 (default) or with the ML model (use_ml = true).
# The match against the species expected from the sample prefix is the same for both scanners.

workflow wf_species_detection_per_sample {

  input {
    File   read1
    File   read2
    String sample_id

    Boolean use_ml = false        # false = kraken2 (default), true = ML model
    File?   ml_model              # optional: replaces the model built into the image (strain_model.pkl)
    Float?  ml_threshold          # empty = threshold stored in the model
    String  docker = "ghcr.io/bioinfomoh/species_detection:2.0"
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
        ml_model = ml_model,          # empty = model built into the image
        sample_id = sample_id,
        threshold = ml_threshold,
        docker = docker
    }
  }

  # Same match logic (sample_id prefix -> expected species) for kraken2 and for the ML model
  call task_match.species_match {
    input:
      sample_id = sample_id,
      detected_genus = select_first([kraken2.genus_detected, ml_species_predict.genus_detected]),
      detected_species = select_first([kraken2.species_detected, ml_species_predict.species_detected])
  }

  output {
    String sp_detection_version          = select_first([kraken2.version, ml_species_predict.version])
    File   sp_detection_report           = select_first([kraken2.report, ml_species_predict.report])
    String sp_detection_species_detected = select_first([kraken2.species_detected, ml_species_predict.species_detected])
    String sp_detection_top_taxa         = select_first([kraken2.top_taxa, ml_species_predict.top_taxa])
    String sp_detection_expected_species = species_match.expected_species
    String sp_detection_match            = species_match.match      # v | ~ | xxx | ??
    String sp_detection_scanner          = if use_ml then "ML" else "kraken2"
  }
}
