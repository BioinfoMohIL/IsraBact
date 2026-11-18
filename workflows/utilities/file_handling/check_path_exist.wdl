version 1.0

# Version "une ligne = un job" pour PETITES tables (quelques centaines/milliers de lignes).
# Lancer sur la table (root entity = ta table, ex. seq_data), avec:
#   read1 = this.read1        read2 = this.read2   (ou laisser vide: déduit de read1)
# et dans l'onglet Outputs:   check_path_exist        = this.check_path_exist          (yes / no)
#                             check_path_exist_report = this.check_path_exist_report   (détail)
# -> Terra écrit directement ces colonnes dans chaque ligne sélectionnée.
# Décoche aussi "Use call caching" à l'envoi (et `volatile: true` ci-dessous l'évite déjà côté WDL).
# Pour de très grandes tables (100 000 lignes), utilise plutôt check_path_exist.wdl (par lots).

workflow check_path_exist_row {
  input {
    String  read1
    String? read2
    String  docker = "google/cloud-sdk:slim"
  }

  call check {
    input: read1 = read1, read2 = read2, docker = docker
  }

  output {
    String check_path_exist        = check.exists   # yes / no
    String check_path_exist_report = check.report   # détail
  }
}

task check {
  input {
    String  read1
    String? read2
    String  docker
  }

  command <<<
    set +e
    r1="~{read1}"
    r2="~{default="" read2}"
    if [ -z "$r2" ]; then
      r2=$(echo "$r1" | sed -e 's/_R1_001\.f/_R2_001.f/' -e 's/_R1\.f/_R2.f/' -e 's/_1\.f/_2.f/')
      [ "$r2" = "$r1" ] && r2=""
    fi
    problems=""
    if [ -z "$r1" ]; then problems="R1_EMPTY"; elif ! gsutil -q stat "$r1"; then problems="R1_MISSING"; fi
    if [ -z "$r2" ]; then problems="${problems:+$problems;}R2_NOT_DERIVABLE"; elif ! gsutil -q stat "$r2"; then problems="${problems:+$problems;}R2_MISSING"; fi
    printf '%s' "${problems:-OK}" > report.txt
    if [ -z "$problems" ]; then printf 'yes' > exists.txt; else printf 'no' > exists.txt; fi
    cat exists.txt; echo; cat report.txt
  >>>

  output {
    String exists = read_string("exists.txt")   # yes / no
    String report = read_string("report.txt")   # OK / R1_MISSING / R2_MISSING / ...
  }

  # Le résultat dépend de l'état du bucket (pas seulement des inputs): ne jamais le servir depuis le cache.
  meta {
    volatile: true
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}
