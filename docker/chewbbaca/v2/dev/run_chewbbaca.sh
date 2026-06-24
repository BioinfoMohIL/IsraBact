#!/bin/bash
set -euo pipefail

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --i) INPUT_ASSEMBLIES="$2"; shift ;;
        --o) OUTPUT_DIR="$2"; shift ;;
        --sample_prefix) SAMPLE_PREFIX="$2"; shift ;;
        --call_dir) CALL_DIR="$2"; shift ;;
        --visualization_dir) VISUALIZATION_DIR="$2"; shift ;;
        *) printf "[ERROR] Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Default paths
WD="/app"
DATASET_DIR="${WD}/datasets"
SCHEMA_PATH="${DATASET_DIR}/${SAMPLE_PREFIX}/${SAMPLE_PREFIX}_scheme.zip"
PRODIGAL_PATH="${DATASET_DIR}/${SAMPLE_PREFIX}/${SAMPLE_PREFIX}_prodigal.zip"
CALL_RESULTS="${CALL_DIR}/*"
SCHEMA_ADAPTED_DIR="${OUTPUT_DIR}/schema_adapted"

echo -e  '\n===== CGMLST Chewbbaca Analysis =====\n'

# 1. Prep External Schema
echo "Unzipping Schema..."
mkdir -p schema_unzipped
unzip "$SCHEMA_PATH" -d schema_unzipped

echo "Unzipping Prodigal..."
mkdir -p prodigal_unzipped
unzip "$PRODIGAL_PATH" -d prodigal_unzipped

PRODIGAL_PATH="${DATASET_DIR}/${SAMPLE_PREFIX}/${SAMPLE_PREFIX}_prodigal.zip"

echo -e  "\n[1] Adapt your external schema"
echo Running: chewBBACA.py PrepExternalSchema -i schema_unzipped -o ${SCHEMA_ADAPTED_DIR} --ptf prodigal_unzipped --cpu 20
chewBBACA.py PrepExternalSchema -i schema_unzipped -o "${SCHEMA_ADAPTED_DIR}" --ptf prodigal_unzipped/*.trn --cpu 20

# 2. Allele Calling
echo -e  "\n[2] Alleles Calling"
echo -e  Running: chewBBACA.py AlleleCall -i ${INPUT_ASSEMBLIES} -g ${SCHEMA_ADAPTED_DIR} -o ${CALL_DIR} --cpu 15
chewBBACA.py AlleleCall -i "${INPUT_ASSEMBLIES}" -g "${SCHEMA_ADAPTED_DIR}" -o "${CALL_DIR}" --cpu 15

# 3. Clean up results_alleles.tsv
echo -e  "\n[3] Cleaning sample name for better tree visualization"

# This will expand to a real path like: results/call/blablabla/results_alleles.tsv
results_file=$(find "$CALL_DIR" -name "results_alleles.tsv" | head -n 1)

awk -F'\t' 'BEGIN{OFS="\t"} {
    name = $1
    sub(/\.[^.]+$/, "", name)              # Remove extension (e.g., .fa, .fasta, .fna)
    split(name, parts, /[_-]/)             # Split on "_" or "-"
    $1 = parts[1]                          # Keep only the first part
    print
}' "$results_file" > ./results_alleles.tsv

# 4. ExtractCgMLST
echo -e  "\n[4] Extraction for data visualization"
echo -e  Running: chewBBACA.py ExtractCgMLST -i ./results_alleles.tsv -o ${VISUALIZATION_DIR}
chewBBACA.py ExtractCgMLST -i ./results_alleles.tsv -o "${VISUALIZATION_DIR}"

echo -e "\n✅ Workflow complete. Output ready for $SAMPLE_PREFIX.\n"
