version 1.0

import "../../../tasks/species_detection/functional_profiling/task_microbedecoder.wdl" as tasks

workflow wf_microbedecoder {
  input {
    String samplename

    # --- Proteins: .faa OR GenBank (.gbff/.gbk), the format is detected automatically ---
    File input_proteins

    # --- KO annotation: "kofamscan" (default) or "eggnog" ---
    String ko_annotator = "kofamscan"
    File? kofam_profiles_tar     # profiles.tar.gz  (if kofamscan)
    File? kofam_ko_list          # ko_list.gz       (if kofamscan)
    File? eggnog_db              # eggnog.db        (if eggnog)
    File? eggnog_dmnd            # eggnog_proteins.dmnd
    File? eggnog_taxa_db         # eggnog.taxa.db

    # --- Metabolic networks (all optional, nothing set = no network is run) ---
    Array[String]? reference_networks        # names of networks from the built-in database (one run per network)
    Array[File]? reference_network_files     # custom network CSV files (one run per file)
    # Shared by all networks: leave empty to use the defaults of each network
    Array[String] substrates = []
    Array[String] products = []
    Array[String] unbalanced_intermediates = []
    Boolean all_subunits = true

    # --- Machine learning (optional: unset or empty = no ML) ---
    # e.g. ["Fermentation (type of metabolism)", "Anaerobe (oxygen tolerance)"]
    Array[String]? model_names

    # --- Taxonomy (optional): set species_name to enable it ---
    String? species_name         # e.g. gambit_predicted_taxon
    String taxonomy_system = "NCBI"
    Array[String] traits_to_predict = ["Type of metabolism (FAPROTAX)"]
    Boolean ignore_NA = true
    Boolean match_all_ranks = false
    Boolean ignore_species = false
  }

  meta {
    description: "Prediction of bacterial metabolic traits (metabolic networks, machine learning, taxonomy) from annotated proteins (Bakta/Prokka). Based on Microbe Decoder (https://github.com/thackmann/MicrobeDecoder.git), run in headless mode."
  }

  # Empty arrays used as fallbacks (select_first) for the optional array inputs
  Array[String] no_strings = []
  Array[File] no_files = []

  call tasks.prepare_proteins {
    input:
      input_proteins = input_proteins,
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

  # Networks from the built-in database (by name)
  scatter (net in select_first([reference_networks, no_strings])) {
    call tasks.microbedecoder_networks as networks_db {
      input:
        gene_functions_file = ko_table,
        reference_network = net,
        substrates = substrates,
        products = products,
        unbalanced_intermediates = unbalanced_intermediates,
        all_subunits = all_subunits,
        samplename = samplename + "_" + sub(net, "[^A-Za-z0-9]+", "_")
    }
  }

  # Custom networks (CSV files)
  scatter (f in select_first([reference_network_files, no_files])) {
    call tasks.microbedecoder_networks as networks_custom {
      input:
        gene_functions_file = ko_table,
        reference_network_file = f,
        substrates = substrates,
        products = products,
        unbalanced_intermediates = unbalanced_intermediates,
        all_subunits = all_subunits,
        samplename = samplename + "_" + sub(basename(f, ".csv"), "[^A-Za-z0-9]+", "_")
    }
  }

  # Optional ML: only runs if at least one model is requested
  if (length(select_first([model_names, no_strings])) > 0) {
    call tasks.microbedecoder_ml {
      input:
        gene_functions_file = ko_table,
        model_names = select_first([model_names, no_strings]),
        samplename = samplename
    }
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
        samplename = samplename
    }
  }

  output {
    File gene_functions_ko = ko_table
    Array[File] networks_fluxes = flatten([networks_db.fluxes_tsv, networks_custom.fluxes_tsv])
    Array[File] networks_rds = flatten([networks_db.result_rds, networks_custom.result_rds])
    Array[File] networks_reactions = flatten([flatten(networks_db.extra_results), flatten(networks_custom.extra_results)])
    File? ml_predictions = microbedecoder_ml.predictions_tsv
    File? ml_rds = microbedecoder_ml.result_rds
    File? taxonomy_predictions = microbedecoder_taxonomy.predictions_tsv
    File? taxonomy_rds = microbedecoder_taxonomy.result_rds
  }
}
