version 1.0

# tasks/task_amr_visualization.wdl
# All pipeline tasks grouped together:
#   - parse_amr   : resistance (AMRFinder/ResFinder) + optional virulence (Abricate)
#                   -> gene x sample matrix
#   - annotate    : Prokka annotation of an assembly  -> GFF3   (AMRViz-like phylo)
#   - pangenome   : Roary pangenome from GFFs          -> core gene alignment
#   - build_tree  : IQ-TREE2 from a core alignment     -> Newick tree
#   - plot_heatmap: phylogeny + resistome/virulome heatmap (PDF/PNG/SVG/HTML)


# ---------------------------------------------------------------------------
task parse_amr {
  input {
    Array[File] resistance_reports               # MANDATORY (AMRFinder or ResFinder)
    Array[File] virulence_reports = []           # OPTIONAL (Abricate vfdb)
    Float min_identity = 90.0
    Float min_coverage = 80.0
    String mode = "presence"
    String out_prefix = "amr"

    String lang = "py"                  # "py" (default) or "r"
    String docker
    String docker_r = "bioinfomoh/amr-heatmap-r:latest"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb = 20
  }

  command <<<
    set -euo pipefail
    ~{if lang == "r" then "parse_amr.R" else "parse_amr.py"} \
      --resistance ~{sep=" " resistance_reports} \
      ~{if length(virulence_reports) > 0 then "--virulence" else ""} ~{sep=" " virulence_reports} \
      --min-identity ~{min_identity} \
      --min-coverage ~{min_coverage} \
      --mode ~{mode} \
      --out ~{out_prefix}_matrix.tsv
  >>>

  output {
    File matrix = "~{out_prefix}_matrix.tsv"
    File genes  = "~{out_prefix}_matrix.tsv.genes.tsv"
  }
  runtime {
    docker: if lang == "r" then docker_r else docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Prokka annotation. The output GFF is named after the (normalized) sample id,
# because Roary derives the sample label from the GFF filename - so this name
# must match the sample ids in the matrix for the tree/heatmap to line up.
task annotate {
  input {
    File assembly
    String sample = sub(basename(assembly), "(_filtered_contigs|_contigs)?\\.(fasta|fa|fna)$", "")
    String docker = "staphb/prokka:1.14.6"
    Int cpu = 2
    Int memory_gb = 4
    Int disk_gb = 20
  }

  command <<<
    set -euo pipefail
    prokka --cpus ~{cpu} --outdir prokka_out --prefix "~{sample}" --force ~{assembly}
    cp prokka_out/"~{sample}".gff ~{sample}.gff
  >>>

  output {
    File gff = "~{sample}.gff"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# Roary pangenome -> core gene alignment (with -e --mafft).
task pangenome {
  input {
    Array[File] gffs
    Int blastp_identity = 95          # roary -i
    String docker = "staphb/roary:3.13.0"
    Int cpu = 4
    Int memory_gb = 16
    Int disk_gb = 50
  }

  command <<<
    set -euo pipefail
    roary -e --mafft -p ~{cpu} -i ~{blastp_identity} -f roary_out ~{sep=" " gffs}
    cp roary_out/core_gene_alignment.aln core_gene_alignment.aln
    cp roary_out/gene_presence_absence.csv gene_presence_absence.csv
  >>>

  output {
    File core_alignment = "core_gene_alignment.aln"
    File gene_presence_absence = "gene_presence_absence.csv"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
# IQ-TREE2 on an EXISTING alignment (Roary core alignment, or a user-supplied
# aligned core_fasta). No MAFFT here - the input is assumed already aligned.
task build_tree {
  input {
    File alignment
    String model = "GTR+G"
    Int bootstrap = 1000
    String out_prefix = "core"
    String docker = "staphb/iqtree2:2.1.2"
    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb = 30
  }

  command <<<
    set -euo pipefail
    iqtree2 -s ~{alignment} -m ~{model} -B ~{bootstrap} -T ~{cpu} --prefix ~{out_prefix}
    cp ~{out_prefix}.treefile ~{out_prefix}.nwk
  >>>

  output {
    File tree       = "~{out_prefix}.nwk"
    File iqtree_log = "~{out_prefix}.iqtree"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


# ---------------------------------------------------------------------------
task plot_heatmap {
  input {
    File matrix
    File genes
    File? tree
    File? metadata
    Array[String] annotate_cols = []
    String title = "Resistome / virulome"
    String out_prefix = "amr_heatmap"

    String lang = "py"                  # "py" (default) or "r"
    String docker
    String docker_r = "bioinfomoh/amr-heatmap-r:latest"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb = 20
  }

  command <<<
    set -euo pipefail
    ~{if lang == "r" then "plot_heatmap.R" else "plot_heatmap.py"} \
      --matrix ~{matrix} \
      --genes ~{genes} \
      ~{if defined(tree) then "--tree " + tree else ""} \
      ~{if defined(metadata) then "--metadata " + metadata else ""} \
      ~{if length(annotate_cols) > 0 then "--annotate" else ""} ~{sep=" " annotate_cols} \
      --title "~{title}" \
      --out ~{out_prefix}.pdf \
      --png ~{out_prefix}.png \
      --svg ~{out_prefix}.svg \
      ~{if lang == "py" then "--html " + out_prefix + ".html" else ""}
  >>>

  output {
    File pdf  = "~{out_prefix}.pdf"
    File png  = "~{out_prefix}.png"
    File svg  = "~{out_prefix}.svg"
    Array[File] html = glob("~{out_prefix}.html")   # py: 1 file; r: empty
  }
  runtime {
    docker: if lang == "r" then docker_r else docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}
