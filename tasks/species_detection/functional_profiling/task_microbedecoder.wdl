version 1.0

# Microbe Decoder (tjhackmann/microbedecoder) en mode headless, apres TheiaProk (Bakta ou Prokka).
# L'image est une app R Shiny : on charge son code avec initialize_session()
# puis on appelle les memes fonctions que l'API Plumber du repo.

# ---------------------------------------------------------------------------
# 1) Proteines : .faa direct, ou extraction depuis un GenBank (.gbff Bakta / .gbk Prokka)
# ---------------------------------------------------------------------------
task prepare_proteins {
  input {
    File? proteins_faa
    File? genbank_file
    String samplename
    String ko_annotator   # verifie ici : "kofamscan" ou "eggnog"

    String docker = "ubuntu:22.04"
    Int cpu = 1
    Int memory_gb = 2
    Int disk_size_gb = 10
  }

  command <<<
    set -euo pipefail

    case "~{ko_annotator}" in
      kofamscan|eggnog) ;;
      *) echo "ko_annotator doit etre 'kofamscan' ou 'eggnog' (recu: ~{ko_annotator})" >&2; exit 1 ;;
    esac

    FAA="~{default="" proteins_faa}"
    GBK="~{default="" genbank_file}"

    if [ -n "$FAA" ]; then
      cp "$FAA" "~{samplename}_proteins.faa"
    elif [ -n "$GBK" ]; then
      awk '
        /^     [^ ]/ { incds = ($1 == "CDS"); intrans = 0; if (incds) locus = ""; next }
        incds && /\/locus_tag=/ { match($0, /"[^"]*"/); locus = substr($0, RSTART + 1, RLENGTH - 2); next }
        incds && /\/translation="/ {
          s = $0; sub(/.*\/translation="/, "", s); gsub(/ /, "", s)
          if (s ~ /"$/) { gsub(/"/, "", s); print ">" locus; print s }
          else { buf = s; intrans = 1 }
          next
        }
        intrans { s = $0; gsub(/ /, "", s)
          if (s ~ /"$/) { gsub(/"/, "", s); print ">" locus; print buf s; intrans = 0 }
          else buf = buf s }
      ' "$GBK" > "~{samplename}_proteins.faa"
    else
      echo "Fournir proteins_faa (.faa de Bakta/Prokka) ou genbank_file (.gbff/.gbk)" >&2
      exit 1
    fi

    [ -s "~{samplename}_proteins.faa" ] || { echo "Aucune proteine trouvee" >&2; exit 1; }
    echo "Proteines: $(grep -c '>' "~{samplename}_proteins.faa")"
  >>>

  output {
    File proteins = "~{samplename}_proteins.faa"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} SSD"
    preemptible: 1
    maxRetries: 1
  }
}

# ---------------------------------------------------------------------------
# 2a) KO avec KofamScan (defaut) : base KOfam = profiles.tar.gz + ko_list(.gz)
#     telechargeables sur ftp://ftp.genome.jp/pub/db/kofam/
# ---------------------------------------------------------------------------
task kofamscan_ko {
  input {
    File proteins_faa
    String samplename
    File kofam_profiles_tar     # profiles.tar.gz
    File kofam_ko_list          # ko_list.gz (ou ko_list decompresse)

    String docker = "quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2"   # a verifier avec docker pull
    Int cpu = 8
    Int memory_gb = 16
    Int disk_size_gb = 50
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail

    mkdir -p kofam
    tar -xzf "~{kofam_profiles_tar}" -C kofam
    case "~{kofam_ko_list}" in
      *.gz) gunzip -c "~{kofam_ko_list}" > kofam/ko_list ;;
      *)    cp "~{kofam_ko_list}" kofam/ko_list ;;
    esac

    HAL="$(find kofam -name prokaryote.hal | head -n 1)"
    [ -n "$HAL" ] || { echo "prokaryote.hal introuvable dans profiles.tar.gz" >&2; exit 1; }

    exec_annotation \
      -p "$HAL" \
      -k kofam/ko_list \
      --cpu ~{cpu} \
      -f mapper \
      -o kofam_mapper.txt \
      "~{proteins_faa}"

    # mapper : "gene<TAB>KO" (les genes sans KO n'ont qu'une colonne)
    awk -F'\t' -v g="~{samplename}" '$2 ~ /^K[0-9][0-9][0-9][0-9][0-9]$/ { print g "\t" $2 }' kofam_mapper.txt | sort -u > ko_rows.tsv

    [ -s ko_rows.tsv ] || { echo "Aucun KO trouve" >&2; exit 1; }
    { printf 'Genome\tKO\n'; cat ko_rows.tsv; } > "~{samplename}_gene_functions.tsv"
  >>>

  output {
    File gene_functions_tsv = "~{samplename}_gene_functions.tsv"
    File kofam_raw = "kofam_mapper.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} SSD"
    preemptible: preemptible
    maxRetries: 1
  }
}

# ---------------------------------------------------------------------------
# 2b) KO avec eggNOG-mapper (option) : sur les proteines (pas de prodigal)
# ---------------------------------------------------------------------------
task eggnog_ko {
  input {
    File proteins_faa
    String samplename
    File eggnog_db              # eggnog.db
    File eggnog_dmnd            # eggnog_proteins.dmnd
    File eggnog_taxa_db         # eggnog.taxa.db

    String docker = "quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_size_gb = 200
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail

    mkdir -p eggnog_data out
    ln -s "~{eggnog_db}" eggnog_data/eggnog.db
    ln -s "~{eggnog_dmnd}" eggnog_data/eggnog_proteins.dmnd
    ln -s "~{eggnog_taxa_db}" eggnog_data/eggnog.taxa.db

    emapper.py \
      -i "~{proteins_faa}" \
      --itype proteins \
      -m diamond \
      --data_dir eggnog_data \
      --output "~{samplename}" --output_dir out \
      --cpu ~{cpu} --override

    awk -F'\t' -v g="~{samplename}" '
      /^#query/ { for (i = 1; i <= NF; i++) if ($i == "KEGG_ko") c = i; next }
      /^#/ { next }
      c && $c != "-" && $c != "" {
        n = split($c, a, ",")
        for (i = 1; i <= n; i++) {
          k = a[i]; sub(/^ko:/, "", k)
          if (k ~ /^K[0-9][0-9][0-9][0-9][0-9]$/) print g "\t" k
        }
      }' "out/~{samplename}.emapper.annotations" | sort -u > ko_rows.tsv

    [ -s ko_rows.tsv ] || { echo "Aucun KO trouve" >&2; exit 1; }
    { printf 'Genome\tKO\n'; cat ko_rows.tsv; } > "~{samplename}_gene_functions.tsv"
  >>>

  output {
    File gene_functions_tsv = "~{samplename}_gene_functions.tsv"
    File emapper_annotations = "out/~{samplename}.emapper.annotations"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} SSD"
    preemptible: preemptible
    maxRetries: 1
  }
}

# ---------------------------------------------------------------------------
# 3) Microbe Decoder - taxonomie (optionnelle) : a partir d'un nom d'espece (GAMBIT)
# ---------------------------------------------------------------------------
task microbedecoder_taxonomy {
  input {
    String species_name                 # ex: gambit_predicted_taxon = "Salmonella enterica"
    Array[String] traits_to_predict = ["Type of metabolism (FAPROTAX)"]
    String system_taxonomy = "NCBI"     # LPSN, GTDB, NCBI ou Bergey
    Boolean ignore_NA = true
    Boolean match_all_ranks = false
    Boolean ignore_species = false
    String query_string = ""
    String samplename = "microbedecoder"

    String docker = "tjhackmann/microbedecoder:latest"
    Int cpu = 2
    Int memory_gb = 8
    Int disk_size_gb = 30
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    OUT="$(pwd)/results"
    mkdir -p "$OUT"

    # Mini-tableau au format DADA2 : seuls Genus et Species (epithete) sont remplis
    SPECIES="~{species_name}"
    GENUS="$(echo "$SPECIES" | awk '{print $1}')"
    EPITHET="$(echo "$SPECIES" | awk '{print $2}')"
    case "$EPITHET" in sp*|"") EPITHET="" ;; esac
    TAXA="$(pwd)/~{samplename}_query_taxa.csv"
    printf 'Phylum,Class,Order,Family,Genus,Species\n,,,,%s,%s\n' "$GENUS" "$EPITHET" > "$TAXA"

    cat > run_taxonomy.R <<'RSCRIPT'
    args <- commandArgs(trailingOnly = TRUE)
    taxa_file <- args[1]; traits_file <- args[2]; system_taxonomy <- args[3]
    ignore_NA <- as.logical(args[4]); match_all_ranks <- as.logical(args[5])
    ignore_species <- as.logical(args[6]); query_string <- args[7]
    out_dir <- args[8]; prefix <- args[9]

    setwd("/srv/shiny-server")
    Sys.setlocale("LC_ALL", "C")
    source("functions/sessionFunctions.R")
    initialize_session(app_dir = getwd(), verbose = FALSE)

    write_table <- function(x, path) {
      df <- tryCatch(as.data.frame(x), error = function(e) NULL)
      if (is.null(df)) writeLines(capture.output(print(x)), path)
      else write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    traits <- trimws(readLines(traits_file, warn = FALSE))
    traits <- traits[nzchar(traits)]

    query <- validate_and_read_file(file_path = taxa_file)
    query[] <- lapply(query, function(x) { x <- trimws(as.character(x)); x[x == ""] <- NA; x })
    query_taxa <- process_uploaded_taxonomy(query)
    if (is.null(query_taxa)) stop("Format de taxonomie non reconnu.")

    result <- compute_taxonomy_predictions(
      data              = load_database(),
      query_taxa        = query_taxa,
      query_string      = if (nzchar(query_string)) query_string else NULL,
      traits_to_predict = traits,
      ignore_NA         = ignore_NA,
      match_all_ranks   = match_all_ranks,
      ignore_species    = ignore_species,
      system_taxonomy   = system_taxonomy
    )

    saveRDS(result, file.path(out_dir, paste0(prefix, "_taxonomy.rds")))
    write_table(result$probabilities, file.path(out_dir, paste0(prefix, "_taxonomy_predictions.tsv")))
    RSCRIPT

    Rscript run_taxonomy.R \
      "$TAXA" \
      "~{write_lines(traits_to_predict)}" \
      "~{system_taxonomy}" \
      "~{ignore_NA}" "~{match_all_ranks}" "~{ignore_species}" \
      "~{query_string}" \
      "$OUT" "~{samplename}"
  >>>

  output {
    File predictions_tsv = "results/~{samplename}_taxonomy_predictions.tsv"
    File result_rds = "results/~{samplename}_taxonomy.rds"
    File query_taxa_csv = "~{samplename}_query_taxa.csv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} SSD"
    preemptible: preemptible
    maxRetries: 1
  }
}

# ---------------------------------------------------------------------------
# 4) Microbe Decoder - reseaux metaboliques
# ---------------------------------------------------------------------------
task microbedecoder_networks {
  input {
    File? gene_functions_file
    Array[String] selected_organisms = []

    String reference_network = "Fermentation of glucose"
    # Si vides, les valeurs par defaut du reseau choisi sont utilisees
    Array[String] substrates = []
    Array[String] products = []
    Array[String] unbalanced_intermediates = []
    Boolean all_subunits = true
    String samplename = "microbedecoder"

    String docker = "tjhackmann/microbedecoder:latest"
    Int cpu = 2
    Int memory_gb = 16
    Int disk_size_gb = 30
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    OUT="$(pwd)/results"
    mkdir -p "$OUT"

    cat > run_networks.R <<'RSCRIPT'
    args <- commandArgs(trailingOnly = TRUE)
    gf_file <- args[1]; organisms_file <- args[2]; reference_network <- args[3]
    substrates_file <- args[4]; products_file <- args[5]; unbalanced_file <- args[6]
    all_subunits <- as.logical(args[7]); out_dir <- args[8]; prefix <- args[9]

    setwd("/srv/shiny-server")
    Sys.setlocale("LC_ALL", "C")
    source("functions/sessionFunctions.R")
    initialize_session(app_dir = getwd(), verbose = FALSE)

    write_table <- function(x, path) {
      df <- tryCatch(as.data.frame(x), error = function(e) NULL)
      if (is.null(df)) writeLines(capture.output(print(x)), path)
      else write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    read_lines_clean <- function(f) { x <- trimws(readLines(f, warn = FALSE)); x[nzchar(x)] }
    pick <- function(f, col) {
      x <- read_lines_clean(f)
      if (length(x) > 0) x else get_metabolite_selections(selected_network = reference_network, metabolite_col = col)
    }

    organisms <- read_lines_clean(organisms_file)
    if (nzchar(gf_file)) {
      gene_functions <- process_uploaded_gene_functions(validate_and_read_file(file_path = gf_file))
      if (is.null(gene_functions)) stop("Format de fonctions genetiques non reconnu.")
    } else if (length(organisms) > 0) {
      gene_functions <- get_gene_functions_from_database(organisms)
    } else {
      stop("Fournir gene_functions_file ou selected_organisms.")
    }

    result <- compute_network_predictions(
      reference_network        = get_reference_network_from_database(reference_network),
      gene_functions           = gene_functions,
      substrates               = pick(substrates_file, "default_substrates"),
      products                 = pick(products_file, "default_products"),
      unbalanced_intermediates = pick(unbalanced_file, "default_unbalanced_intermediates"),
      all_subunits             = all_subunits
    )

    saveRDS(result, file.path(out_dir, paste0(prefix, "_networks.rds")))
    write_table(result$fluxes, file.path(out_dir, paste0(prefix, "_networks_fluxes.tsv")))
    if (is.data.frame(result$reaction_table)) {
      write_table(result$reaction_table, file.path(out_dir, paste0(prefix, "_networks_reactions.tsv")))
    }
    RSCRIPT

    Rscript run_networks.R \
      "~{default="" gene_functions_file}" \
      "~{write_lines(selected_organisms)}" \
      "~{reference_network}" \
      "~{write_lines(substrates)}" \
      "~{write_lines(products)}" \
      "~{write_lines(unbalanced_intermediates)}" \
      "~{all_subunits}" \
      "$OUT" "~{samplename}"
  >>>

  output {
    File fluxes_tsv = "results/~{samplename}_networks_fluxes.tsv"
    File result_rds = "results/~{samplename}_networks.rds"
    Array[File] extra_results = glob("results/*_networks_reactions.tsv")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} SSD"
    preemptible: preemptible
    maxRetries: 1
  }
}

# ---------------------------------------------------------------------------
# 5) Microbe Decoder - machine learning
# ---------------------------------------------------------------------------
task microbedecoder_ml {
  input {
    File? gene_functions_file
    Array[String] selected_organisms = []
    # Noms exacts des modeles, ex: "Fermentation (type of metabolism)", "Anaerobe (oxygen tolerance)"
    Array[String] model_names = ["Fermentation (type of metabolism)", "Methanogenesis (type of metabolism)"]
    String samplename = "microbedecoder"

    String docker = "tjhackmann/microbedecoder:latest"
    Int cpu = 2
    Int memory_gb = 16
    Int disk_size_gb = 30
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    OUT="$(pwd)/results"
    mkdir -p "$OUT"

    cat > run_ml.R <<'RSCRIPT'
    args <- commandArgs(trailingOnly = TRUE)
    gf_file <- args[1]; organisms_file <- args[2]; models_file <- args[3]
    out_dir <- args[4]; prefix <- args[5]

    setwd("/srv/shiny-server")
    Sys.setlocale("LC_ALL", "C")
    source("functions/sessionFunctions.R")
    initialize_session(app_dir = getwd(), verbose = FALSE)

    write_table <- function(x, path) {
      df <- tryCatch(as.data.frame(x), error = function(e) NULL)
      if (is.null(df)) writeLines(capture.output(print(x)), path)
      else write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    read_lines_clean <- function(f) { x <- trimws(readLines(f, warn = FALSE)); x[nzchar(x)] }

    organisms   <- read_lines_clean(organisms_file)
    model_names <- read_lines_clean(models_file)
    if (nzchar(gf_file)) {
      gene_functions <- process_uploaded_gene_functions(validate_and_read_file(file_path = gf_file))
      if (is.null(gene_functions)) stop("Format de fonctions genetiques non reconnu.")
    } else if (length(organisms) > 0) {
      gene_functions <- get_gene_functions_from_database(organisms)
    } else {
      stop("Fournir gene_functions_file ou selected_organisms.")
    }

    model_paths <- get_model_paths_from_database(model_names, model_path_config = model_path_config)

    result <- compute_ml_predictions(
      df = gene_functions, model_names = model_names, model_paths = model_paths,
      response = NULL, predictors = NULL, seed = NULL, ntree = NULL, maxnodes = NULL,
      positive_class_weight = NULL, training_split = NULL, keep_models = FALSE
    )

    saveRDS(result, file.path(out_dir, paste0(prefix, "_ml.rds")))
    write_table(result$probabilities, file.path(out_dir, paste0(prefix, "_ml_predictions.tsv")))
    RSCRIPT

    Rscript run_ml.R \
      "~{default="" gene_functions_file}" \
      "~{write_lines(selected_organisms)}" \
      "~{write_lines(model_names)}" \
      "$OUT" "~{samplename}"
  >>>

  output {
    File predictions_tsv = "results/~{samplename}_ml_predictions.tsv"
    File result_rds = "results/~{samplename}_ml.rds"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size_gb} SSD"
    preemptible: preemptible
    maxRetries: 1
  }
}
