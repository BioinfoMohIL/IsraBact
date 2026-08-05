version 1.0

## task_build_gene_matrix.wdl
##
## Prend en entree une liste de rapports TSV (AbriCate OU AMRFinderPlus, un fichier
## par souche) et construit 2 matrices "souche x gene" :
##   - une matrice texte YES/NO (presence/absence)
##   - une matrice binaire 0/1
## Chaque matrice est produite en TSV et en XLSX.
##
## La cle (Key = nom de la souche) est deduite du nom du fichier d'entree
## (tout ce qui precede le premier "_"), ex:
##   V0022_amrfinder_all.tsv   -> V0022
##   CA116948_abricate_hits.tsv -> CA116948

task build_gene_matrix {
  input {
    Array[File] tsv_reports     # un fichier de rapport par souche
    String gene_column          # nom exact de la colonne "gene" dans le TSV
                                 #   AbriCate      -> "GENE"
                                 #   AMRFinderPlus -> "Element symbol"
    String source_name          # prefixe utilise pour nommer les fichiers de sortie
                                 #   ex: "abricate", "amrfinder"

    Int mem_gb = 4
    Int cpu = 1
    Int disk_gb = 20
    String docker = "python:3.11-slim"
  }

  command <<<
    set -euo pipefail

    pip install --quiet --no-cache-dir pandas openpyxl

    python3 <<'PYEOF'
import os
import pandas as pd

files = """~{sep='\n' tsv_reports}""".strip().split("\n")
gene_col = "~{gene_column}"
source_name = "~{source_name}"

presence = {}   # Key -> set(genes)
all_genes = set()

for f in files:
    f = f.strip()
    if not f:
        continue
    basename = os.path.basename(f)
    key = basename.split("_")[0]

    df = pd.read_csv(f, sep="\t", dtype=str, comment=None)
    df.columns = [c.strip().lstrip("#") for c in df.columns]

    if gene_col not in df.columns:
        raise ValueError(
            f"Colonne '{gene_col}' introuvable dans {basename}. "
            f"Colonnes disponibles: {list(df.columns)}"
        )

    genes = set(df[gene_col].dropna().astype(str).str.strip())
    genes.discard("")

    presence.setdefault(key, set())
    presence[key] |= genes
    all_genes |= genes

keys_sorted = sorted(presence.keys())
genes_sorted = sorted(all_genes)

yesno_rows = []
binary_rows = []
for k in keys_sorted:
    yesno_row = {"Key": k}
    binary_row = {"Key": k}
    for g in genes_sorted:
        present = g in presence[k]
        yesno_row[g] = "YES" if present else "NO"
        binary_row[g] = 1 if present else 0
    yesno_rows.append(yesno_row)
    binary_rows.append(binary_row)

cols = ["Key"] + genes_sorted
yesno_df = pd.DataFrame(yesno_rows, columns=cols)
binary_df = pd.DataFrame(binary_rows, columns=cols)

yesno_df.to_csv(f"{source_name}_matrix_yesno.tsv", sep="\t", index=False)
binary_df.to_csv(f"{source_name}_matrix_binary.tsv", sep="\t", index=False)

yesno_df.to_excel(f"{source_name}_matrix_yesno.xlsx", index=False)
binary_df.to_excel(f"{source_name}_matrix_binary.xlsx", index=False)
PYEOF
  >>>

  output {
    File matrix_yesno_tsv   = "~{source_name}_matrix_yesno.tsv"
    File matrix_binary_tsv  = "~{source_name}_matrix_binary.tsv"
    File matrix_yesno_xlsx  = "~{source_name}_matrix_yesno.xlsx"
    File matrix_binary_xlsx = "~{source_name}_matrix_binary.xlsx"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: cpu
    disks: "local-disk " + disk_gb + " HDD"
  }
}
