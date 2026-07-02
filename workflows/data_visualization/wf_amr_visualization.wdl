version 1.0

# workflows/wf_amr_visualization.wdl
# AMRViz-like pipeline: AMR reports -> matrix, and (optionally) assemblies ->
# Prokka -> Roary -> IQ-TREE -> tree, then a phylogeny + heatmap figure.

import "../../tasks/data_visualization/task_amr_visualization.wdl" as tasks

workflow wf_amr_visualization {
  input {
    # --- AMR inputs ---
    Array[File] resistance_reports               # MANDATORY (AMRFinder or ResFinder)
    Array[File] virulence_reports = []           # OPTIONAL (Abricate vfdb)
    Float min_identity = 90.0
    Float min_coverage = 80.0
    String mode = "presence"

    # --- phylogeny: pick ONE of three options (checked in this order) ---
    File? tree                 # 1) precomputed Newick -> used as-is
    File? core_fasta           # 2) an existing core ALIGNMENT -> IQ-TREE only
    Array[File] assemblies = [] # 3) raw assemblies -> Prokka -> Roary -> IQ-TREE
    String iqtree_model = "GTR+G"
    Int bootstrap = 1000
    Int roary_identity = 95

    # --- annotation / rendering ---
    File? metadata
    Array[String] annotate_cols = []
    String title = "Resistome / virulome"
    String collection = "amr"

    # --- language for parse + plot: "py" (default) or "r" ---
    String lang = "r"

    # --- docker images ---
    String docker                                # python image (parse + plot)
    String docker_r = "bioinfomoh/amr-heatmap-r:latest"   # R image (ggtree)
    String prokka_docker = "staphb/prokka:1.14.6"
    String roary_docker  = "staphb/roary:3.13.0"
    String iqtree_docker = "staphb/iqtree2:2.1.2"
  }

  # 1) gene x sample matrix
  call tasks.parse_amr {
    input:
      resistance_reports = resistance_reports,
      virulence_reports = virulence_reports,
      min_identity = min_identity,
      min_coverage = min_coverage,
      mode = mode,
      out_prefix = collection,
      lang = lang,
      docker = docker,
      docker_r = docker_r
  }

  # 2) AMRViz-like phylogeny from assemblies (only if no tree/alignment given)
  if (!defined(tree) && !defined(core_fasta) && length(assemblies) > 0) {
    scatter (asm in assemblies) {
      call tasks.annotate as annotate_asm {
        input: assembly = asm, docker = prokka_docker
      }
    }
    call tasks.pangenome {
      input: gffs = annotate_asm.gff, blastp_identity = roary_identity, docker = roary_docker
    }
  }

  # alignment to feed IQ-TREE: a user-supplied core alignment, or Roary's output
  File? alignment = if defined(core_fasta) then core_fasta else pangenome.core_alignment

  if (!defined(tree) && defined(alignment)) {
    call tasks.build_tree {
      input:
        alignment = select_first([alignment]),
        model = iqtree_model,
        bootstrap = bootstrap,
        out_prefix = collection,
        docker = iqtree_docker
    }
  }

  # effective tree: provided > computed > none (then heatmap sorts alphabetically)
  File? effective_tree = if defined(tree) then tree else build_tree.tree

  # 3) figure
  call tasks.plot_heatmap {
    input:
      matrix = parse_amr.matrix,
      genes = parse_amr.genes,
      tree = effective_tree,
      metadata = metadata,
      annotate_cols = annotate_cols,
      title = title,
      out_prefix = collection + "_heatmap",
      lang = lang,
      docker = docker,
      docker_r = docker_r
  }

  output {
    File matrix = parse_amr.matrix
    File genes = parse_amr.genes
    File? tree_used = effective_tree
    File? gene_presence_absence = pangenome.gene_presence_absence
    File heatmap_pdf = plot_heatmap.pdf
    File heatmap_png = plot_heatmap.png
    File heatmap_svg = plot_heatmap.svg
    Array[File] heatmap_html = plot_heatmap.html   # py: 1 HTML; r: empty
  }
}
