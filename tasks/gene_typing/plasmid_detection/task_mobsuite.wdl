version 1.0

## Task: MOB-suite (mob_recon)
## Reconstructs and types plasmids from a draft assembly (contigs FASTA).
## Reference: https://github.com/phac-nml/mob-suite
## Docker: staphb/mob-suite:3.1.9
##
## Every parameter that has a literal default value in mob_recon.py (v3.1.9)
## is exposed here with that same default. Parameters with no simple default
## (filter files, database overrides, boolean flags) are optional and are
## only added to the command line if provided.

task mob_recon {
  input {
    File   assembly_fasta       # -i/--infile : assembled contigs (FASTA) for the sample
    String samplename           # -s/--sample_id

    # --- Parameters with an official default (mob_recon.py, v3.1.9) ---
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

    # --- Boolean flags (default = absent = False in mob_recon) ---
    Boolean unicycler_contigs = false   # -u : circularity via unicycler header
    Boolean run_overhang      = false   # -c : circularity detection via assembly overhang
    Boolean keep_tmp          = false   # -k : keep the temporary directory
    Boolean debug              = false  # --debug

    # --- Optionals with no simple default (omitted from the command if unset) ---
    File?   filter_db                  # -b : fasta of sequences to mask
    String? genome_filter_db_prefix    # -g : mash sketch / blastdb prefix for closed genomes
    String? prefix                     # -p : prefix for result files
    String? database_directory         # -d : override for the MOB-suite database directory
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

    # mob_recon can exit non-zero on some messy assemblies (no plasmid, contigs
    # too short, etc.). We capture the return code instead of letting a single
    # failed sample break the whole workflow.
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

    # Normalize expected output filenames even if mob_recon failed before
    # generating everything, so the WDL outputs stay valid (empty files).
    touch ~{out_dir}/contig_report.txt
    touch ~{out_dir}/mobtyper_results.txt
    touch ~{out_dir}/mge.report.txt
    touch ~{out_dir}/chromosome.fasta
    # biomarkers.blast.txt is NOT touched here: mob_recon only writes it when
    # plasmid biomarkers are actually detected (otherwise it never creates the
    # file at all). We leave its absence as a true "0 biomarkers found" signal
    # instead of masking it with an artificial empty file.

    # Number of plasmids detected = number of plasmid_*.fasta files
    PLASMID_COUNT=$(ls ~{out_dir}/plasmid_*.fasta 2>/dev/null | wc -l || true)
    echo "$PLASMID_COUNT" > PLASMID_COUNT

    # Bundle all plasmid fasta files into a single tar archive, even if empty
    tar -czf ~{samplename}_plasmids.tar.gz -C ~{out_dir} $(cd ~{out_dir} && ls plasmid_*.fasta 2>/dev/null) 2>/dev/null || tar -czf ~{samplename}_plasmids.tar.gz --files-from=/dev/null

    # Rename the tabular text reports to .tsv: a .txt file's text/plain MIME
    # type gets opened inline by most browsers when downloaded from Terra,
    # while .tsv (text/tab-separated-values) generally forces a real download.
    mv ~{out_dir}/contig_report.txt ~{out_dir}/contig_report.tsv
    mv ~{out_dir}/mobtyper_results.txt ~{out_dir}/mobtyper_results.tsv
    mv ~{out_dir}/mge.report.txt ~{out_dir}/mge_report.tsv
    if [ -f ~{out_dir}/biomarkers.blast.txt ]; then
      mv ~{out_dir}/biomarkers.blast.txt ~{out_dir}/biomarkers_blast.tsv
    fi

    if [ "$MOB_EXIT" -eq 0 ]; then
      echo "SUCCESS" > MOB_RECON_STATUS
    elif [ "$MOB_EXIT" -eq 255 ]; then
        echo "FAILED" > MOB_RECON_STATUS
        echo "mob_recon returned 255 - likely transient database initialization/download error" >&2

        exit 255
    else
        echo "FAILED" > MOB_RECON_STATUS
        echo "mob_recon returned a non-zero exit code (~{samplename}): $MOB_EXIT" >&2
    fi
  >>>

  output {
    String mob_suite_version = read_string("VERSION")
    String mob_recon_status  = read_string("MOB_RECON_STATUS")  # "SUCCESS" or "FAILED"
    Int    plasmid_count     = read_int("PLASMID_COUNT")

    File   contig_report       = "~{out_dir}/contig_report.tsv"
    File   mobtyper_report     = "~{out_dir}/mobtyper_results.tsv"
    File   mge_report          = "~{out_dir}/mge_report.tsv"
    File?  biomarker_report    = "~{out_dir}/biomarkers_blast.tsv"  # null if mob_recon found nothing
    File   chromosome_fasta    = "~{out_dir}/chromosome.fasta"
    File   plasmids_tarball    = "~{samplename}_plasmids.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 1
    maxRetries: 2
  }

  meta {
    description: "MOB-suite (mob_recon) reconstructs and types plasmid sequences from draft or complete bacterial genome assemblies - (https://github.com/phac-nml/mob-suite)."
    author: "David Maimoun"
  }
}
