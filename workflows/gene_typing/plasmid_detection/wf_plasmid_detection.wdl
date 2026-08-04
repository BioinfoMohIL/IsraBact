version 1.0

import "task_mobsuite.wdl" as mobsuite_t

## Workflow: Plasmid_Detection_PerSample_MOH
## Detection et typage de plasmides via MOB-suite (mob_recon) a partir d'un
## assemblage de contigs. Concu pour Terra/GCP Batch, single-click, un sample
## a la fois (comme wf_species_detection_bs).
## Reference outil: https://github.com/phac-nml/mob-suite

workflow wf_plasmid_detection {
  input {
    File   assembly_fasta
    String samplename

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

    Boolean unicycler_contigs = false
    Boolean run_overhang      = false
    Boolean keep_tmp          = false
    Boolean debug              = false

    File?   filter_db
    String? genome_filter_db_prefix
    String? prefix
    String? database_directory
    File?   plasmid_db
    File?   plasmid_mash_db
    File?   plasmid_meta
    String? plasmid_db_type
    File?   plasmid_replicons
    File?   repetitive_mask
    File?   plasmid_mob
    File?   plasmid_mpf
    File?   plasmid_orit

    Int    cpu       = 4
    Int    memory    = 8
    Int    disk_size = 50
    String docker    = "staphb/mob-suite:3.1.9"
  }

  call mobsuite_t.mob_recon {
    input:
      assembly_fasta                  = assembly_fasta,
      samplename                      = samplename,
      mash_genome_neighbor_threshold  = mash_genome_neighbor_threshold,
      max_contig_size                 = max_contig_size,
      max_plasmid_size                = max_plasmid_size,
      min_rep_evalue                  = min_rep_evalue,
      min_mob_evalue                  = min_mob_evalue,
      min_con_evalue                  = min_con_evalue,
      min_rpp_evalue                  = min_rpp_evalue,
      min_length                      = min_length,
      min_rep_ident                   = min_rep_ident,
      min_mob_ident                   = min_mob_ident,
      min_con_ident                   = min_con_ident,
      min_rpp_ident                   = min_rpp_ident,
      min_rep_cov                     = min_rep_cov,
      min_mob_cov                     = min_mob_cov,
      min_con_cov                     = min_con_cov,
      min_rpp_cov                     = min_rpp_cov,
      min_overlap                     = min_overlap,
      primary_cluster_dist            = primary_cluster_dist,
      secondary_cluster_dist          = secondary_cluster_dist,
      unicycler_contigs                = unicycler_contigs,
      run_overhang                     = run_overhang,
      keep_tmp                         = keep_tmp,
      debug                            = debug,
      filter_db                        = filter_db,
      genome_filter_db_prefix          = genome_filter_db_prefix,
      prefix                           = prefix,
      database_directory               = database_directory,
      plasmid_db                       = plasmid_db,
      plasmid_mash_db                  = plasmid_mash_db,
      plasmid_meta                     = plasmid_meta,
      plasmid_db_type                  = plasmid_db_type,
      plasmid_replicons                = plasmid_replicons,
      repetitive_mask                  = repetitive_mask,
      plasmid_mob                      = plasmid_mob,
      plasmid_mpf                      = plasmid_mpf,
      plasmid_orit                     = plasmid_orit,
      cpu                               = cpu,
      memory                            = memory,
      disk_size                         = disk_size,
      docker                            = docker
  }

  output {
    String mob_suite_version   = mob_recon.mob_suite_version
    Int    mob_recon_exit_code = mob_recon.mob_recon_exit_code
    Int    plasmid_count       = mob_recon.plasmid_count

    File   contig_report       = mob_recon.contig_report
    File   mobtyper_report     = mob_recon.mobtyper_report
    File   mge_report          = mob_recon.mge_report
    File?  biomarker_report    = mob_recon.biomarker_report
    File   chromosome_fasta    = mob_recon.chromosome_fasta
    File   plasmids_tarball    = mob_recon.plasmids_tarball
  }

  meta {
    description: "Detection et typage de plasmides via MOB-suite (mob_recon), un sample a la fois."
  }
}
