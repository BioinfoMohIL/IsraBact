#!/usr/bin/env bash
# Lance le pipeline complet sur les donnees d'exemple (8 isolats K. pneumoniae).
# Usage : bash examples/run_example.sh
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
IN="$HERE/inputs"
OUT="$HERE/outputs"
mkdir -p "$OUT"

echo "[1/2] Parse Abricate + VFDB + ResFinder -> matrice"
python "$ROOT/scripts/parse_amr.py" \
  --resistance "$IN"/resistance/*.tsv \
  --virulence "$IN"/virulence/*.tsv \
  --min-identity 90 --min-coverage 80 --mode presence \
  --out "$OUT/KpMDR_matrix.tsv"

echo "[2/2] Figure phylogenie + heatmap"
python "$ROOT/scripts/plot_heatmap.py" \
  --matrix "$OUT/KpMDR_matrix.tsv" \
  --genes  "$OUT/KpMDR_matrix.tsv.genes.tsv" \
  --tree   "$IN/tree.nwk" \
  --metadata "$IN/metadata.tsv" --annotate ST species \
  --title "MDR K. pneumoniae" \
  --out "$OUT/KpMDR_heatmap.pdf" --png "$OUT/KpMDR_heatmap.png" \
  --svg "$OUT/KpMDR_heatmap.svg" --html "$OUT/KpMDR_heatmap.html"

echo "OK -> $OUT/KpMDR_heatmap.pdf"
