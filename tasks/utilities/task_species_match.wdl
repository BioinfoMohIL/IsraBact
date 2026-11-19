version 1.0

# Compare l'espèce détectée (par kraken2 OU par le ML) à l'espèce attendue d'après le préfixe
# du nom d'échantillon. Logique extraite telle quelle de task_kraken2.
#
#   match = "v"    : conforme à l'attendu
#           "xxx"  : différent
#           "??"   : indéterminable (préfixe inconnu, ou rien détecté)
#
# Attendu sans espace (ex. "Salmonella")              -> comparé au GENRE détecté
# Attendu avec espace (ex. "Listeria monocytogenes")  -> comparé à l'ESPÈCE détectée
# Préfixe: le plus long d'abord (VIB avant V, F-EC avant EC), sur le nom en majuscules.

task species_match {
  input {
    String sample_id
    String detected_genus      # "" si rien n'a été détecté (ou si le ML répond UNKNOWN)
    String detected_species    # nom complet de l'espèce détectée
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
    elif [[ "$expected" == *" "* ]]; then
        # species-level expectation
        if [ "$species" == "$expected" ]; then match="v"; else match="xxx"; fi
    else
        # genus-level expectation
        if [ "$genus" == "$expected" ]; then match="v"; else match="xxx"; fi
    fi

    echo "$match" > match.txt
    echo "${expected:-NA}" > expected.txt
  >>>

  output {
    String match            = read_string("match.txt")
    String expected_species = read_string("expected.txt")   # "NA" si le préfixe est inconnu
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}
