version 1.0

## task_build_gene_matrix.wdl
##
## Takes as input a list of TSV reports (AbriCate OR AMRFinderPlus - amr/stress/
## virulence, one file per strain) and builds:
##   - a "strain x gene" YES/NO matrix (presence/absence)
##   - a "strain x gene" binary 0/1 matrix
##   - a file concatenating all raw rows from the input reports
##     (with a "Key" column added as the first column)
## Each matrix is produced as TSV and XLSX.
##
## The key (Key = strain name) is derived from the input file name
## (everything before the first "_"), e.g.:
##   V0022_amrfinder_amr_report.tsv -> V0022
##   CA116948_abricate_hits.tsv     -> CA116948

task build_gene_matrix {
  input {
    Array[File] tsv_reports     # one report file per strain
    String gene_column          # exact name of the "gene" column in the TSV
                                 #   AbriCate      -> "GENE"
                                 #   AMRFinderPlus -> "Element symbol"
    String source_name          # prefix used to name the output files
                                 #   e.g.: "abricate", "amrfinder_amr", "amrfinder_stress",
                                 #       "amrfinder_virulence"

  }

  command <<<
    set -euo pipefail

    pip install --quiet --no-cache-dir pandas openpyxl

    cat > build_gene_matrix.py <<'PYEOF'
import os
import pandas as pd

files = """~{sep='\n' tsv_reports}""".strip().split("\n")
gene_col = "~{gene_column}"
source_name = "~{source_name}"

presence = {}   # Key -> set(genes)
all_genes = set()
concat_frames = []

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

    df.insert(0, "Key", key)
    concat_frames.append(df)

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

# yesno_df.to_csv(f"{source_name}_matrix_yesno.tsv", sep="\t", index=False)
# binary_df.to_csv(f"{source_name}_matrix_binary.tsv", sep="\t", index=False)

yesno_df.to_excel(f"{source_name}_matrix_yesno.xlsx", index=False)
binary_df.to_excel(f"{source_name}_matrix_binary.xlsx", index=False)

# file concatenating all raw rows from the input reports
concat_df = pd.concat(concat_frames, ignore_index=True, sort=False)
concat_df.to_csv(f"{source_name}_concatenated.tsv", sep="\t", index=False)
PYEOF

    # Run through textwrap.dedent so the script still works even if the
    # WDL command block above ever gets re-indented by an editor/formatter -
    # a uniform leading-whitespace shift would otherwise break Python syntax.
    python3 - <<'BOOT'
import textwrap
with open("build_gene_matrix.py") as fh:
    src = fh.read()
exec(compile(textwrap.dedent(src), "build_gene_matrix.py", "exec"))
BOOT
  >>>

  output {
    # File matrix_yesno_tsv    = "~{source_name}_matrix_yesno.tsv"
    # File matrix_binary_tsv   = "~{source_name}_matrix_binary.tsv"
    File matrix_yesno_xlsx   = "~{source_name}_matrix_yesno.xlsx"
    File matrix_binary_xlsx  = "~{source_name}_matrix_binary.xlsx"
    File concatenated_tsv    = "~{source_name}_concatenated.tsv"
  }

    runtime {
    docker: "python:3.11-slim"
    memory: "4 GB"
    cpu: 4
    disks: "local-disk 20 HDD"
  }
}
