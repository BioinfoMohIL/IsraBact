version 1.0


task mlst {
  input {
    File assembly_fasta
    String taxon
  }

  command <<<
    #!/usr/bin/env bash
  
    pubmlst_db="https://rest.pubmlst.org/db"
    pasteur_db="https://bigsdb.pasteur.fr/api/db"

    taxon=~{taxon}
    taxon_chosen="$(echo "~{taxon}" | tr '[:upper:]' '[:lower:]')"

    is_ec_or_sh=0

    case "$taxon_chosen" in
      "bordetella pertussis")      db="${pasteur_db}/pubmlst_bordetella_seqdef/schemes/3/sequence" ;;
      "listeria monocytogenes")    db="${pasteur_db}/pubmlst_listeria_seqdef/schemes/2/sequence" ;;
      "haemophilus influenzae")    db="${pubmlst_db}/pubmlst_hinfluenzae_seqdef/schemes/1/sequence" ;;
      "neisseria meningitidis")    db="${pubmlst_db}/pubmlst_neisseria_seqdef/schemes/1/sequence" ;;
      "staphylococcus aureus")     db="${pubmlst_db}/pubmlst_saureus_seqdef/schemes/1/sequence" ;;
      "streptococcus pyogenes")    db="${pubmlst_db}/pubmlst_spyogenes_seqdef/schemes/1/sequence" ;;
      "streptococcus agalactiae")  db="${pubmlst_db}/pubmlst_sagalactiae_seqdef/schemes/1/sequence" ;;
      "streptococcus pneumoniae")  db="${pubmlst_db}/pubmlst_spneumoniae_seqdef/schemes/1/sequence" ;;
      "streptococcus intermedius") db="${pubmlst_db}/pubmlst_spseudintermedius_seqdef/schemes/1/sequence" ;;
      "streptococcus dentisani")   db="${pubmlst_db}/pubmlst_oralstrep_seqdef/schemes/1/sequence" ;;
      "vibrio cholerae")           db="${pubmlst_db}/pubmlst_vcholerae_seqdef/schemes/1/sequence" ;;
      "vibrio vulnificus")         db="${pubmlst_db}/pubmlst_vvulnificus_seqdef/schemes/1/sequence" ;;
      "escherichia coli"|"shigella sonnei")  
        db="${pubmlst_db}/pubmlst_mlst_seqdef/schemes/4/sequence"
        is_ec_or_sh=1
        ;;
      "shigella boydii"|"shigella flexneri")  
        db="${pubmlst_db}/pubmlst_mlst_seqdef/schemes/4/sequence"
        is_ec_or_sh=1
        ;;
      "campylobacter jejuni"|"campylobacter coli") db="${pubmlst_db}/pubmlst_campylobacter_seqdef/schemes/1/sequence" ;;
      "salmonella"|"salmonella enterica")          db="${pubmlst_db}/pubmlst_salmonella_seqdef/schemes/2/sequence" ;;
      *)
        
      echo "Unknown taxon name: '$taxon_chosen'" >&2
      exit 2
      ;;
    esac

    if [[ "$db" =~ (/|^)(pubmlst_[^/]+_seqdef)(/|$) ]]; then 
      echo "${BASH_REMATCH[2]}" > db_used.txt
    fi

    (echo -n '{"base64":true,"sequence":"'; base64 "~{assembly_fasta}"; echo '"}') | \
      curl -s -H "Content-Type: application/json" \
      -X POST "${db}" -d @- > mlst.json

    if [[ ! -s mlst.json ]]; then
        echo "❌ Error: mlst.json not found or empty"
        exit 1
    fi

    jq -r 'if .fields.ST then "ST\(.fields.ST)" else empty end' mlst.json > st.txt

    # Remove Escherichia_coli_Achtman_MLST_ form all alleles
    prefix_to_remove=""
    if [[ "$is_ec_or_sh" -eq 1 ]]; then
        prefix_to_remove="Escherichia_coli_Achtman_MLST_"
    fi

    jq -r --arg prefix "$prefix_to_remove" '
      (.exact_matches // {}) 
      | to_entries 
      | map(
          "\(.key)(\(.value[0].allele_id))" 
          | sub($prefix;"")
        )
      | sort
      | join(",")
    ' mlst.json > alleles.txt

  >>>

  output {
    File result     = "mlst.json"
    String db_used  = read_string("db_used.txt")
    String alleles  = read_string('alleles.txt')
    String st       = read_string('st.txt')
  }

  runtime {
    docker: "bioinfomoh/linux-tools:1"
    cpu: 1
    memory: "1G"
  }
}
