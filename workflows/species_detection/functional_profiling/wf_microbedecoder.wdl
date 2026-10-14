version 1.0

import "../../../tasks/species_detection/functional_profiling/task_microbedecoder.wdl" as tasks

# Microbe Decoder apres TheiaProk (annotation Bakta ou Prokka).
#   Metabolique (toujours) : proteines -> KO (KofamScan par defaut, ou eggNOG) -> reseaux + ML
#   Taxonomie (optionnelle) : nom d'espece (GAMBIT) -> Microbe Decoder
#
# Entree proteines : .faa de Bakta/Prokka (proteins_faa), sinon un GenBank (.gbff/.gbk).
# Le TSV de Bakta ne contient pas les sequences : il ne peut pas servir d'entree.

workflow wf_microbedecoder {
  input {
    String samplename
    String docker = "tjhackmann/microbedecoder:latest"

    # --- Proteines (au moins un des deux) ---
    File? proteins_faa
    File? genbank_file

    # --- Annotation KO : "kofamscan" (defaut) ou "eggnog" ---
    String ko_annotator = "kofamscan"
    File? kofam_profiles_tar     # profiles.tar.gz  (si kofamscan)
    File? kofam_ko_list          # ko_list.gz       (si kofamscan)
    File? eggnog_db              # eggnog.db        (si eggnog)
    File? eggnog_dmnd            # eggnog_proteins.dmnd
    File? eggnog_taxa_db         # eggnog.taxa.db

    # --- Taxonomie (optionnelle) : renseigner species_name pour l'activer ---
    String? species_name         # ex: gambit_predicted_taxon
    String taxonomy_system = "NCBI"
    Array[String] traits_to_predict = ["Type of metabolism (FAPROTAX)"]
    Boolean ignore_NA = true
    Boolean match_all_ranks = false
    Boolean ignore_species = false

    # --- Reseaux metaboliques ---
    String reference_network = "Fermentation of glucose"
    Array[String] substrates = []
    Array[String] products = []
    Array[String] unbalanced_intermediates = []
    Boolean all_subunits = true

    # --- Machine learning ---
    Array[String] model_names = ["Fermentation (type of metabolism)", "Methanogenesis (type of metabolism)"]
  }

  call tasks.prepare_proteins {
    input:
      proteins_faa = proteins_faa,
      genbank_file = genbank_file,
      samplename = samplename,
      ko_annotator = ko_annotator
  }

  if (ko_annotator == "kofamscan") {
    call tasks.kofamscan_ko {
      input:
        proteins_faa = prepare_proteins.proteins,
        samplename = samplename,
        kofam_profiles_tar = select_first([kofam_profiles_tar]),
        kofam_ko_list = select_first([kofam_ko_list])
    }
  }

  if (ko_annotator == "eggnog") {
    call tasks.eggnog_ko {
      input:
        proteins_faa = prepare_proteins.proteins,
        samplename = samplename,
        eggnog_db = select_first([eggnog_db]),
        eggnog_dmnd = select_first([eggnog_dmnd]),
        eggnog_taxa_db = select_first([eggnog_taxa_db])
    }
  }

  File ko_table = select_first([kofamscan_ko.gene_functions_tsv, eggnog_ko.gene_functions_tsv])

  call tasks.microbedecoder_networks {
    input:
      gene_functions_file = ko_table,
      reference_network = reference_network,
      substrates = substrates,
      products = products,
      unbalanced_intermediates = unbalanced_intermediates,
      all_subunits = all_subunits,
      samplename = samplename,
      docker = docker
  }

  call tasks.microbedecoder_ml {
    input:
      gene_functions_file = ko_table,
      model_names = model_names,
      samplename = samplename,
      docker = docker
  }

  if (defined(species_name)) {
    call tasks.microbedecoder_taxonomy {
      input:
        species_name = select_first([species_name]),
        traits_to_predict = traits_to_predict,
        system_taxonomy = taxonomy_system,
        ignore_NA = ignore_NA,
        match_all_ranks = match_all_ranks,
        ignore_species = ignore_species,
        samplename = samplename,
        docker = docker
    }
  }

  output {
    File gene_functions_ko = ko_table
    File networks_fluxes = microbedecoder_networks.fluxes_tsv
    File networks_rds = microbedecoder_networks.result_rds
    File ml_predictions = microbedecoder_ml.predictions_tsv
    File ml_rds = microbedecoder_ml.result_rds
    File? taxonomy_predictions = microbedecoder_taxonomy.predictions_tsv
    File? taxonomy_rds = microbedecoder_taxonomy.result_rds
  }
}
