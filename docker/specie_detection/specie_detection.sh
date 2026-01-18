#!/bin/bash

# Set default values
cpu=1
read1=""
read2=""
output_report="kraken_report.txt"
specie_detected="sd.txt"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --read1)
      read1="$2"
      shift 2
      ;;
    --read2)
      read2="$2"
      shift 2
      ;;
    --cpu)
      cpu="$2"
      shift 2
      ;;
    --output_report)
      output_report="$2"
      shift 2
      ;;
    --specie_detected)
      specie_detected="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

# Check if required parameters are provided
if [ -z "$read1" ]; then
  echo "Error: --read1 is required"
  exit 1
fi

# Initialize variables for kraken2 options
mode=""
compressed=""

# Check if paired mode should be used
if ! [ -z "$read2" ]; then
    echo "Reads are paired..."
    mode="--paired"
fi

# Determine if reads are compressed
if [[ "$read1" == *.gz ]]; then
    echo "Reads are compressed..."
    compressed="--gzip-compressed"
fi

# Run Kraken2
echo "Running Kraken2..."
kraken2 $mode $compressed --threads "$cpu" --use-names --db /app/db/kraken_db \
    --report "$output_report" \
    --paired "$read1" "$read2" \
    --output -

# Extract taxonomy information
awk -F'\t' '$4 == "S" {gsub(/^[ \t]+/, "", $6); print $6; exit}' "$output_report" > "$specie_detected"