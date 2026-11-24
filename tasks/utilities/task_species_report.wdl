version 1.0

# One report for a whole run (a set of samples) and a one-line summary of it, instead of
# returning long per-sample arrays.
#
#   species_report    : arrays (one value per sample)  ->  a single CSV
#                         sample_id,detected_species,match,expected_species,scanner,top_taxa
#                       match: v = as expected | ~ = same genus, other species (approx.) | xxx = different | ?? = unknown
#   summarize_species : that CSV  ->  "ALL MATCHED!" or e.g. "MATCH: 90,APPROX: 2,NOMATCH: 3,UNKNOWN: 1"
#                       (adapted from SummarizeSpecies: "+" is now "v", and "~" has its own count, APPROX)

task species_report {
  input {
    Array[String] sample_ids
    Array[String] species_detected    # same order as sample_ids
    Array[String] matches
    Array[String] expected_species
    Array[String] top_taxa
    String        scanner             # "kraken2" or "ML"
    String        output_name = "species_detection_report.csv"
    String        docker = "ubuntu:22.04"
  }

  command <<<
    set -euo pipefail

    paste ~{write_lines(sample_ids)} ~{write_lines(species_detected)} ~{write_lines(matches)} \
          ~{write_lines(expected_species)} ~{write_lines(top_taxa)} > rows.tsv

    # match stays in column 3 (like "sample,detected,match" in the previous workflow).
    # Only fields that contain a comma or a quote are quoted (top_taxa: "Genus:98.7, Species:35.0").
    echo "sample_id,detected_species,match,expected_species,scanner,top_taxa" > "~{output_name}"
    awk -F'\t' -v scanner="~{scanner}" '
      function q(s) { if (s ~ /[",]/) { gsub(/"/, "\"\"", s); return "\"" s "\"" } return s }
      { print q($1) "," q($2) "," q($3) "," q($4) "," scanner "," q($5) }
    ' rows.tsv >> "~{output_name}"
  >>>

  output {
    File report = "~{output_name}"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}

task summarize_species {
  input {
    File   report_csv                 # output of species_report (the "match" column is found by its header)
    String docker = "ubuntu:22.04"
  }

  command <<<
    set -euo pipefail

    counts=$(awk -F',' '
      NR == 1 { for (i = 1; i <= NF; i++) if ($i == "match") c = i; next }
      c { if ($c == "v") m++; else if ($c == "~") a++; else if ($c == "xxx") x++; else if ($c == "??") u++ }
      END { print m + 0, a + 0, x + 0, u + 0 }
    ' "~{report_csv}")
    read -r MATCH APPROX NOMATCH UNKNOWN <<< "$counts"

    if [ "$MATCH" -gt 0 ] && [ "$APPROX" -eq 0 ] && [ "$NOMATCH" -eq 0 ] && [ "$UNKNOWN" -eq 0 ]; then
        echo "ALL MATCHED!" > summary.txt
    else
        SUMMARY=""
        add() { if [ "$2" -gt 0 ]; then SUMMARY="${SUMMARY:+$SUMMARY,}$1: $2"; fi; }
        add MATCH "$MATCH"
        add APPROX "$APPROX"
        add NOMATCH "$NOMATCH"
        add UNKNOWN "$UNKNOWN"
        echo "${SUMMARY:-NO SAMPLES}" > summary.txt
    fi
  >>>

  output {
    String summary = read_string("summary.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}
