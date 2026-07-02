#!/usr/bin/env Rscript
# parse_amr.R - R twin of parse_amr.py. Same CLI and outputs.
#   --resistance f1 f2 ...   (mandatory; AMRFinder or ResFinder, auto-detected)
#   --virulence  f1 f2 ...   (optional; Abricate vfdb)
#   --min-identity --min-coverage --mode {presence,identity} --out
# Base R only (no extra packages needed for parsing).

STRIP <- c("_resfinder",".resfinder","_abricate",".abricate","_amrfinder",
  "_amr","_all","_hits","_results_tab","_pheno_table","_phenotype",
  ".tsv",".txt",".tab",".tab.txt","_filtered_contigs","_filtered","_contigs",
  "_assembly","_spades","_skesa","_shovill","_polished",
  ".fasta",".fa",".fna",".contigs",".assembly")

normalize_id <- function(raw) {
  name <- basename(trimws(as.character(raw)))
  repeat {
    changed <- FALSE; low <- tolower(name)
    for (s in STRIP) if (endsWith(low, s)) {
      name <- substr(name, 1, nchar(name) - nchar(s)); changed <- TRUE; low <- tolower(name)
    }
    if (!changed) break
  }
  if (nchar(name) == 0) basename(as.character(raw)) else name
}

detect_format <- function(path) {
  ln <- readLines(path, n = 40, warn = FALSE); ne <- ln[nchar(trimws(ln)) > 0]
  if (length(ne) == 0) return("unknown")
  if (startsWith(tolower(ne[1]), "# resfinder phenotype")) return("resfinder_pheno")
  for (l in ne) {
    low <- tolower(l)
    if (startsWith(low, "#file\t") || startsWith(low, "#file ")) return("abricate")
    if (!startsWith(l, "#")) {
      cols <- tolower(trimws(strsplit(l, "\t")[[1]]))
      if ("element symbol" %in% cols) return("amrfinder")
      if ("resistance gene" %in% cols) return("resfinder_tab")
      if ("gene" %in% cols && "database" %in% cols) return("abricate")
      break
    }
  }
  "unknown"
}

.col <- function(df, nm) { i <- which(tolower(trimws(names(df))) == nm); if (length(i)) df[[i[1]]] else NULL }
.num <- function(x, d = 100) { v <- suppressWarnings(as.numeric(x)); v[is.na(v)] <- d; v }
EMPTY <- function() data.frame(sample = character(), gene = character(),
  identity = numeric(), category = character(), drug_class = character(), stringsAsFactors = FALSE)

parse_amrfinder <- function(path, minid, mincov) {
  df <- read.delim(path, check.names = FALSE, colClasses = "character", quote = "")
  smp <- vapply(.col(df, "name"), normalize_id, ""); gene <- trimws(.col(df, "element symbol"))
  ident <- .num(.col(df, "% identity to reference")); cov <- .num(.col(df, "% coverage of reference"))
  typ <- toupper(trimws(.col(df, "type"))); cls <- .col(df, "class"); if (is.null(cls)) cls <- rep("", nrow(df))
  cmap <- c(AMR = "resistance", VIRULENCE = "virulence", STRESS = "stress")
  cat <- ifelse(typ %in% names(cmap), cmap[typ], "resistance")
  keep <- ident >= minid & cov >= mincov & nchar(gene) > 0
  list(rec = data.frame(sample = smp[keep], gene = gene[keep], identity = ident[keep],
       category = unname(cat[keep]), drug_class = cls[keep], stringsAsFactors = FALSE),
       samples = unique(smp))
}

parse_resfinder_pheno <- function(path, minid, mincov) {
  sample <- normalize_id(path); rows <- list()
  for (line in readLines(path, warn = FALSE)) {
    if (startsWith(line, "# Sample:")) sample <- normalize_id(sub("^# Sample:", "", line))
    if (startsWith(line, "#") || nchar(trimws(line)) == 0) next
    p <- trimws(strsplit(line, "\t")[[1]]); if (length(p) < 5) next
    if (nchar(p[5]) == 0) next
    for (tok in strsplit(p[5], ",")[[1]]) {
      g <- trimws(strsplit(tok, " \\(")[[1]][1])
      if (nchar(g)) rows[[length(rows) + 1]] <- data.frame(sample = sample, gene = g,
        identity = 100, category = "resistance", drug_class = p[2], stringsAsFactors = FALSE)
    }
  }
  list(rec = if (length(rows)) do.call(rbind, rows) else EMPTY(), samples = sample)
}

parse_resfinder_tab <- function(path, minid, mincov) {
  df <- read.delim(path, check.names = FALSE, colClasses = "character", quote = "")
  gene <- trimws(.col(df, "resistance gene")); if (is.null(gene)) gene <- trimws(.col(df, "gene"))
  ident <- .num(.col(df, "identity")); cov <- .num(.col(df, "coverage"))
  cls <- .col(df, "phenotype"); if (is.null(cls)) cls <- rep("", nrow(df))
  sample <- normalize_id(path); keep <- ident >= minid & cov >= mincov & nchar(gene) > 0
  list(rec = data.frame(sample = rep(sample, sum(keep)), gene = gene[keep], identity = ident[keep],
       category = "resistance", drug_class = cls[keep], stringsAsFactors = FALSE), samples = sample)
}

parse_abricate <- function(path, minid, mincov, default_cat = "virulence") {
  df <- read.delim(path, check.names = FALSE, colClasses = "character", quote = "")
  names(df) <- toupper(sub("^#", "", trimws(names(df))))
  if (!"GENE" %in% names(df)) return(list(rec = EMPTY(), samples = character()))
  smp <- if ("FILE" %in% names(df)) vapply(df[["FILE"]], normalize_id, "") else rep(normalize_id(path), nrow(df))
  ident <- .num(df[["%IDENTITY"]]); cov <- .num(df[["%COVERAGE"]])
  db <- tolower(trimws(if ("DATABASE" %in% names(df)) df[["DATABASE"]] else rep("", nrow(df))))
  dbmap <- c(vfdb = "virulence", ecoli_vf = "virulence", resfinder = "resistance",
             card = "resistance", ncbi = "resistance", argannot = "resistance", plasmidfinder = "plasmid")
  cat <- ifelse(db %in% names(dbmap), dbmap[db], default_cat)
  keep <- ident >= minid & cov >= mincov
  list(rec = data.frame(sample = smp[keep], gene = trimws(df[["GENE"]])[keep], identity = ident[keep],
       category = unname(cat[keep]), drug_class = "", stringsAsFactors = FALSE), samples = unique(smp))
}

parse_resistance <- function(p, mi, mc) {
  f <- detect_format(p)
  if (f == "amrfinder") parse_amrfinder(p, mi, mc)
  else if (f == "resfinder_pheno") parse_resfinder_pheno(p, mi, mc)
  else if (f == "resfinder_tab") parse_resfinder_tab(p, mi, mc)
  else if (f == "abricate") parse_abricate(p, mi, mc, "resistance")
  else { message("[warn] unrecognized resistance: ", p); list(rec = EMPTY(), samples = character()) }
}
parse_virulence <- function(p, mi, mc) {
  f <- detect_format(p)
  if (f == "abricate") parse_abricate(p, mi, mc, "virulence")
  else if (f == "amrfinder") { r <- parse_amrfinder(p, mi, mc); r$rec <- r$rec[r$rec$category == "virulence", ]; r }
  else { message("[warn] unrecognized virulence: ", p); list(rec = EMPTY(), samples = character()) }
}

build_matrix <- function(rec, all_samples, mode) {
  samples <- sort(unique(c(rec$sample, all_samples)))
  if (nrow(rec) == 0) return(list(mat = matrix(numeric(0), length(samples), 0, dimnames = list(samples, NULL)),
                                   ann = data.frame(gene = character(), category = character(), drug_class = character())))
  rec <- rec[order(rec$identity), ]
  rec <- rec[!duplicated(rec[, c("sample", "gene")], fromLast = TRUE), ]
  rec$val <- if (mode == "presence") 1 else rec$identity
  genes <- sort(unique(rec$gene))
  m <- matrix(0, length(samples), length(genes), dimnames = list(samples, genes))
  m[cbind(match(rec$sample, samples), match(rec$gene, genes))] <- rec$val
  cat_of <- tapply(rec$category, rec$gene, function(x) names(sort(table(x), decreasing = TRUE))[1])
  dc_of <- tapply(rec$drug_class, rec$gene, function(x) { x <- x[x != ""]; if (length(x)) x[1] else "" })
  ann <- data.frame(gene = genes, category = cat_of[genes], drug_class = dc_of[genes], stringsAsFactors = FALSE)
  list(mat = m, ann = ann)
}

## --- manual CLI parser (handles list-valued flags like --resistance f1 f2) ---
ARGS <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, multiple = FALSE, default = NULL) {
  i <- which(ARGS == name)
  if (!length(i)) return(if (multiple) character() else default)
  vals <- character(); j <- i[1] + 1
  while (j <= length(ARGS) && !startsWith(ARGS[j], "--")) { vals <- c(vals, ARGS[j]); j <- j + 1 }
  if (multiple) vals else if (length(vals)) vals[1] else default
}

res_files <- getflag("--resistance", multiple = TRUE)
vir_files <- getflag("--virulence", multiple = TRUE)
min_identity <- as.numeric(getflag("--min-identity", default = "90"))
min_coverage <- as.numeric(getflag("--min-coverage", default = "80"))
mode <- getflag("--mode", default = "presence")
out <- getflag("--out")
if (!length(res_files)) stop("--resistance is required")

rec <- EMPTY(); alls <- character()
for (p in res_files) { r <- parse_resistance(p, min_identity, min_coverage); rec <- rbind(rec, r$rec); alls <- union(alls, r$samples) }
for (p in vir_files) { r <- parse_virulence(p, min_identity, min_coverage); rec <- rbind(rec, r$rec); alls <- union(alls, r$samples) }
if (length(alls) == 0) stop("No samples parsed: check resistance inputs and thresholds.")

bm <- build_matrix(rec, alls, mode)
out_df <- data.frame(sample = rownames(bm$mat), bm$mat, check.names = FALSE)
write.table(out_df, out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bm$ann, paste0(out, ".genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[ok] %d samples x %d genes -> %s", nrow(bm$mat), ncol(bm$mat), out))
