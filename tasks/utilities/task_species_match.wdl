version 1.0

# Compares the detected species (kraken2 OR the ML model) with the species expected from the
# prefix of the sample name. Shared by both scanners.
#
#   match = "v"    : as expected
#           "~"    : approximately: same genus but a different species (kraken2 is not always precise)
#           "xxx"  : different
#           "??"   : cannot be decided (unknown prefix, nothing detected / ML answered UNKNOWN,
#                  or Shigella detected as E. coli: kraken2 confuses them, not an error of ours)
#
# Expected name without a space (e.g. "Salmonella")             -> compared with the detected GENUS
# Expected name with a space (e.g. "Listeria monocytogenes")    -> compared with the detected SPECIES
#                                                                   (same genus only -> "~")
# Prefix: longest first (VIB before V, F-EC before EC), on the upper-cased sample name.

task species_match {
  input {
    String sample_id
    String detected_genus      # "" if nothing was detected (or the ML model answered UNKNOWN)
    String detected_species    # full name of the detected species
    String docker = "ubuntu:22.04"
  }

  command <<<
    set -euo pipefail

    genus="~{detected_genus}"
    species="~{detected_species}"

    # --- Expected species from sample prefix ---
    declare -A EXPECTED=(
        ["SH"]="Salmonella"
        ["SO"]="Salmonella"
        ["LC"]="Listeria monocytogenes"
        ["LF"]="Listeria monocytogenes"
        ["SG"]="Shigella"
        ["CA"]="Campylobacter"
        ["VIB"]="Vibrio"
        ["V"]="Vibrio"
        ["EC"]="Escherichia coli"
        ["F-EC"]="Escherichia coli"
        ["SA"]="Staphylococcus aureus"
        ["BP"]="Bordetella pertussis"
        ["ST"]="Streptococcus"
        ["SP"]="Streptococcus pneumoniae"
        ["LG"]="Legionella pneumophila"
        ["LW"]="Legionella pneumophila"
        ["CB"]="Corynebacterium diphtheriae"
        ["HI"]="Haemophilus influenzae"
        ["NM"]="Neisseria meningitidis"
        ["M"]="Neisseria meningitidis"
    )

    sid="~{sample_id}"
    sid="${sid^^}"
    expected=""

    # Longest prefix first (so VIB wins over V, F-EC over EC...)
    for k in $(printf '%s\n' "${!EXPECTED[@]}" | awk '{print length($0), $0}' | sort -rn | cut -d' ' -f2); do
        if [[ "$sid" == "$k"* ]]; then
            expected="${EXPECTED[$k]}"
            break
        fi
    done

    # --- Match logic ---
    if [ -z "$expected" ] || [ -z "$genus" ]; then
        match="??"
    elif [ "$expected" == "Shigella" ] && [ "$species" == "Escherichia coli" ]; then
        match="??"         # Shigella is detected as E. coli (known confusion), not an error of ours
    elif [[ "$expected" == *" "* ]]; then
        # species-level expectation
        if [ "$species" == "$expected" ]; then
            match="v"
        elif [ "$genus" == "${expected%% *}" ]; then
            match="~"      # same genus, different species
        else
            match="xxx"
        fi
    else
        # genus-level expectation
        if [ "$genus" == "$expected" ]; then match="v"; else match="xxx"; fi
    fi

    echo "$match" > match.txt
    echo "${expected:-NA}" > expected.txt
  >>>

  output {
    String match            = read_string("match.txt")
    String expected_species = read_string("expected.txt")   # "NA" if the prefix is unknown
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}
