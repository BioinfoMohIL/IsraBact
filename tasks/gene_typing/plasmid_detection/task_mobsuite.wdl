version 1.0

## Task: MOB-suite (mob_recon)
## Reconstruction et typage de plasmides a partir d'un assemblage draft (contigs FASTA).
## Reference: https://github.com/phac-nml/mob-suite
## Docker: staphb/mob-suite:3.1.9
##
## Tous les parametres avec une valeur par defaut litterale dans mob_recon.py
## (v3.1.9) sont exposes ici avec cette meme valeur en dur. Les parametres sans
## defaut simple (fichiers de filtre, overrides de base de donnees, flags
## booleens) sont optionnels et ne sont ajoutes a la commande que si fournis.

task mob_recon {
  input {
    File   assembly_fasta       # -i/--infile : contigs assembles (FASTA) du sample
    String samplename           # -s/--sample_id

    # --- Parametres avec defaut officiel (mob_recon.py, v3.1.9) ---
    Float  mash_genome_neighbor_threshold = 0.002
    Int    max_contig_size                = 450000
    Int    max_plasmid_size               = 450000
    Float  min_rep_evalue                 = 0.00001
    Float  min_mob_evalue                 = 0.00001
    Float  min_con_evalue                 = 0.00001
    Float  min_rpp_evalue                 = 0.00001
    Int    min_length                     = 1000
    Int    min_rep_ident                  = 80
    Int    min_mob_ident                  = 80
    Int    min_con_ident                  = 80
    Int    min_rpp_ident                  = 80
    Int    min_rep_cov                    = 80
    Int    min_mob_cov                    = 80
    Int    min_con_cov                    = 60
    Int    min_rpp_cov                    = 80
    Int    min_overlap                    = 10
    Float  primary_cluster_dist           = 0.06
    Float  secondary_cluster_dist         = 0.025

    # --- Flags booleens (defaut = absent = False chez mob_recon) ---
    Boolean unicycler_contigs = false   # -u : circularite via header unicycler
    Boolean run_overhang      = false   # -c : detection circularite par overhang
    Boolean keep_tmp          = false   # -k : garder le dossier temporaire
    Boolean debug              = false  # --debug

    # --- Optionnels sans defaut simple (non passes si non fournis) ---
    File?   filter_db                  # -b : fasta de sequences a masquer
    String? genome_filter_db_prefix    # -g : prefixe sketch mash/blastdb genomes fermes
    String? prefix                     # -p : prefixe des fichiers de resultats
    String? database_directory         # -d : override du repertoire de bases MOB-suite
    File?   plasmid_db                 # --plasmid_db
    File?   plasmid_mash_db            # --plasmid_mash_db
    File?   plasmid_meta               # -m/--plasmid_meta
    String? plasmid_db_type            # --plasmid_db_type
    File?   plasmid_replicons          # --plasmid_replicons
    File?   repetitive_mask            # --repetitive_mask
    File?   plasmid_mob                # --plasmid_mob
    File?   plasmid_mpf                # --plasmid_mpf
    File?   plasmid_orit               # --plasmid_orit

    Int    cpu       = 4
    Int    memory    = 8
    Int    disk_size = 50
    String docker    = "staphb/mob-suite:3.1.9"
  }

  String out_dir = "~{samplename}_results"

  command <<<
    set -euo pipefail

    mob_recon --version | tee VERSION

    mkdir -p ~{out_dir}

    # mob_recon echoue (exit != 0) sur certains assemblages degueres (pas de plasmide,
    # contigs trop courts, etc.) : on capture le code retour pour ne pas casser
    # tout le workflow sur un seul sample defaillant.
    set +e
    mob_recon \
      --infile ~{assembly_fasta} \
      --outdir ~{out_dir} \
      --sample_id ~{samplename} \
      --num_threads ~{cpu} \
      --mash_genome_neighbor_threshold ~{mash_genome_neighbor_threshold} \
      --max_contig_size ~{max_contig_size} \
      --max_plasmid_size ~{max_plasmid_size} \
      --min_rep_evalue ~{min_rep_evalue} \
      --min_mob_evalue ~{min_mob_evalue} \
      --min_con_evalue ~{min_con_evalue} \
      --min_rpp_evalue ~{min_rpp_evalue} \
      --min_length ~{min_length} \
      --min_rep_ident ~{min_rep_ident} \
      --min_mob_ident ~{min_mob_ident} \
      --min_con_ident ~{min_con_ident} \
      --min_rpp_ident ~{min_rpp_ident} \
      --min_rep_cov ~{min_rep_cov} \
      --min_mob_cov ~{min_mob_cov} \
      --min_con_cov ~{min_con_cov} \
      --min_rpp_cov ~{min_rpp_cov} \
      --min_overlap ~{min_overlap} \
      --primary_cluster_dist ~{primary_cluster_dist} \
      --secondary_cluster_dist ~{secondary_cluster_dist} \
      ~{true="--unicycler_contigs" false="" unicycler_contigs} \
      ~{true="--run_overhang" false="" run_overhang} \
      ~{true="--keep_tmp" false="" keep_tmp} \
      ~{true="--debug" false="" debug} \
      ~{"--filter_db " + filter_db} \
      ~{"--genome_filter_db_prefix " + genome_filter_db_prefix} \
      ~{"--prefix " + prefix} \
      ~{"--database_directory " + database_directory} \
      ~{"--plasmid_db " + plasmid_db} \
      ~{"--plasmid_mash_db " + plasmid_mash_db} \
      ~{"--plasmid_meta " + plasmid_meta} \
      ~{"--plasmid_db_type " + plasmid_db_type} \
      ~{"--plasmid_replicons " + plasmid_replicons} \
      ~{"--repetitive_mask " + repetitive_mask} \
      ~{"--plasmid_mob " + plasmid_mob} \
      ~{"--plasmid_mpf " + plasmid_mpf} \
      ~{"--plasmid_orit " + plasmid_orit} \
      --force
    MOB_EXIT=$?
    set -e

    echo "$MOB_EXIT" > MOB_RECON_EXIT_CODE

    # Normalise les noms de sortie attendus meme si mob_recon a echoue avant
    # de tout generer, pour que les outputs WDL restent valides (fichiers vides).
    touch ~{out_dir}/contig_report.txt
    touch ~{out_dir}/mobtyper_results.txt
    touch ~{out_dir}/mge.report.txt
    touch ~{out_dir}/biomarkers.blast.txt
    touch ~{out_dir}/chromosome.fasta

    # Nombre de plasmides detectes = nombre de fichiers plasmid_*.fasta
    PLASMID_COUNT=$(ls ~{out_dir}/plasmid_*.fasta 2>/dev/null | wc -l || true)
    echo "$PLASMID_COUNT" > PLASMID_COUNT

    # Regroupe tous les fasta de plasmides en une seule archive tar, meme si vide
    tar -czf ~{samplename}_plasmids.tar.gz -C ~{out_dir} $(cd ~{out_dir} && ls plasmid_*.fasta 2>/dev/null) 2>/dev/null || tar -czf ~{samplename}_plasmids.tar.gz --files-from=/dev/null

    if [ "$MOB_EXIT" -ne 0 ]; then
      echo "mob_recon a retourne un code d'erreur (~{samplename}): $MOB_EXIT" >&2
    fi
  >>>

  output {
    String mob_suite_version   = read_string("VERSION")
    Int    mob_recon_exit_code = read_int("MOB_RECON_EXIT_CODE")
    Int    plasmid_count       = read_int("PLASMID_COUNT")

    File   contig_report       = "~{out_dir}/contig_report.txt"
    File   mobtyper_report     = "~{out_dir}/mobtyper_results.txt"
    File   mge_report          = "~{out_dir}/mge.report.txt"
    File   biomarker_report    = "~{out_dir}/biomarkers.blast.txt"
    File   chromosome_fasta    = "~{out_dir}/chromosome.fasta"
    File   plasmids_tarball    = "~{samplename}_plasmids.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 1
    maxRetries: 1
    continueOnReturnCode: true
  }
}
