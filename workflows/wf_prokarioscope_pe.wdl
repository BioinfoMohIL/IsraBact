version 1.0

import "../tasks/task_versioning.wdl" as versioning
import "../tasks/qc/validation/validate_reads.wdl" as validate_reads
import "../tasks/qc/read/task_cg_pipeline.wdl" as cg
import "../tasks/taxon_id/task_gambit.wdl" as gambit_task
import "../tasks/taxon_id/contamination/task_kmerfinder.wdl" as kmerfinder_task
import "../tasks/gene_typing/drug_resistance/task_amrfinderplus.wdl" as amrfinderplus
import "../tasks/gene_typing/drug_resistance/task_resfinder.wdl" as resfinder
import "../tasks/gene_typing/multi/task_ts_mlst.wdl" as ts_mlst_task
import "../tasks/quality_control/advanced_metrics/task_mummer_ani.wdl" as ani_task
import "./utilities/wf_assembly.wdl" as deno_assembly
import "./utilities/wf_read_qc_trim_pe.wdl" as read_QC
import "./utilities/species_typing/wf_pathotype.wdl" as pathotype
import "./utilities/file_handling/wf_concatenate_illumina_lanes.wdl" as concatenate_lanes_workflow


workflow prokarioscope_pe {
    meta {
        description: "De-novo genome assembly, taxonomic ID, and QC of paired-end bacterial NGS data"
        author: "David Maimoun (The Codon Bleu)"
        email: "thecodonbleu@outlook.com"
    }
  
    input {
        String samplename
        String seq_method = "ILLUMINA"
        File read1
        File read2

        # If nextseq lanes
        File? read1_lane2
        File? read1_lane3
        File? read1_lane4
        File? read2_lane2
        File? read2_lane3
        File? read2_lane4
    
    	Int? genome_length

        # Export taxon table parameters
        String? run_id
        String? collection_date
        String? originating_lab
        String? city
        String? county
        String? zip
        File? taxon_tables
        String terra_project = "NA"
        String terra_workspace = "NA"
        
        # read validation parameters
        Boolean is_validate_reads = false 
        Int min_reads = 7472
        Int min_basepairs = 2241820
        Int min_genome_length = 100000
        Int max_genome_length = 18040666
        Int min_coverage = 10
        Int min_proportion = 40
        
        # trimming parameters
        Int trim_min_length = 75
        Int trim_quality_min_score = 20
        Int trim_window_size = 10

        # optional parameters for unicycler (long reads assembly)
        File? illumina_unpaired_fq
        File? long_reads
    
        # module options
        Boolean perform_characterization = true # by default run all characterization steps
        Boolean call_ani = false # by default do not call ANI task, but user has ability to enable this task if working with enteric pathogens or supply their own high-quality reference genome
        Boolean call_kmerfinder = false 
        Boolean call_resfinder = false
        Boolean call_plasmidfinder = true
        Boolean call_abricate = false
        Boolean call_arln_stats = false
        String abricate_db = "vfdb"
        String genome_annotation = "prokka" # options: "prokka" or "bakta"
        String bakta_db = "full" # Default: "light" or "full"
        String? expected_taxon  # allow user to provide organism (e.g. "Clostridioides_difficile") string to amrfinder. Useful when gambit does not predict the correct species    # qc check parameters
        File? qc_check_table

    }

    call versioning.version_capture {
        input:
    }

    if (defined(read1_lane2)) {
        call concatenate_lanes_workflow.concatenate_illumina_lanes {
        input:
            samplename = samplename,
            read1_lane1 = read1,
            read1_lane2 = select_first([read1_lane2]),
            read1_lane3 = read1_lane3,
            read1_lane4 = read1_lane4,
            read2_lane1 = read2,
            read2_lane2 = read2_lane2,
            read2_lane3 = read2_lane3,
            read2_lane4 = read2_lane4
        }
    }


    if (is_validate_reads) {
        call validate_reads.check_reads as raw_check_reads {
        input:
            read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
            read2 = select_first([concatenate_illumina_lanes.read2_concatenated, read2]),
            min_reads = min_reads,
            min_basepairs = min_basepairs,
            min_genome_length = min_genome_length,
            max_genome_length = max_genome_length,
            min_coverage = min_coverage,
            min_proportion = min_proportion,
            expected_genome_length = genome_length
        }
    }

    call read_QC.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
        read2 = select_first([concatenate_illumina_lanes.read2_concatenated, read2]),
        trim_min_length = trim_min_length,
        trim_quality_min_score = trim_quality_min_score,
        trim_window_size = trim_window_size,
        workflow_series = "prokarioscope"
    }

    call deno_assembly.assembly {
        input:
            samplename = samplename,
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            illumina_unpaired_fq = illumina_unpaired_fq,
            long_reads = long_reads
    }

	call cg.cg_pipeline as cg_pipeline_raw {
		input:
			read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
			read2 = select_first([concatenate_illumina_lanes.read2_concatenated, read2]),
			samplename = samplename,
			genome_length = select_first([genome_length, assembly.assembly_length])
	}

	call cg.cg_pipeline as cg_pipeline_clean {
		input:
			read1 = read_QC_trim.read1_clean,
			read2 = read_QC_trim.read2_clean,
			samplename = samplename,
			genome_length = select_first([genome_length, quast.genome_length])
	}

	if (perform_characterization) {
		call gambit_task.gambit {
			input:
			assembly = assembly.assembly_fasta,
			samplename = samplename
		}

		if (call_ani) {
			call ani_task.animummer as ani {
			input:
				assembly = assembly.assembly_fasta,
				samplename = samplename
			}
		}

		if (call_kmerfinder) {
			call kmerfinder_task.kmerfinder_bacteria as kmerfinder {
			input:
				assembly = assembly.assembly_fasta,
				samplename = samplename
			}
		}

		call amrfinderplus.amrfinderplus_nuc as amrfinderplus_task {
			input:
			assembly = assembly.assembly_fasta,
			samplename = samplename,
			organism = select_first([expected_taxon, gambit.gambit_predicted_taxon])
		}

		if (call_resfinder) {
			call resfinder.resfinder as resfinder_task {
			input:
				assembly = assembly.assembly_fasta,
				samplename = samplename,
				organism = select_first([expected_taxon, gambit.gambit_predicted_taxon])
			}
		}

		call ts_mlst_task.ts_mlst {
			input: 
			assembly = assembly.assembly_fasta,
			samplename = samplename,
			taxonomy = select_first([expected_taxon, gambit.gambit_predicted_taxon])
		}
	}

   
    output {
        String prokarioscope_pe_version = version_capture.version
        String prokarioscope_pe_date    = version_capture.date
        
        String seq_platform = seq_method
        
        # If Concatenated Reads
        File? read1_concatenated = concatenate_illumina_lanes.read1_concatenated
    
        # Read QC - outputs
        Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
        Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
        String? fastq_scan_num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
        String? fastq_scan_version = read_QC_trim.fastq_scan_version
        Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
        Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
        String? fastq_scan_num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
        File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
        File? fastq_scan_raw2_json = read_QC_trim.fastq_scan_raw2_json
        File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
        File? fastq_scan_clean2_json = read_QC_trim.fastq_scan_clean2_json
        
        # Read QC - fastqc outputs
        Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
        Int? fastqc_num_reads_raw2 = read_QC_trim.fastqc_raw2
        String? fastqc_num_reads_raw_pairs = read_QC_trim.fastqc_raw_pairs
        Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
        Int? fastqc_num_reads_clean2 = read_QC_trim.fastqc_clean2
        String? fastqc_num_reads_clean_pairs = read_QC_trim.fastqc_clean_pairs
        File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
        File? fastqc_raw2_html = read_QC_trim.fastqc_raw2_html
        File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
        File? fastqc_clean2_html = read_QC_trim.fastqc_clean2_html
        String? fastqc_version = read_QC_trim.fastqc_version
        String? fastqc_docker = read_QC_trim.fastqc_docker
        
        # Read QC - trimmomatic outputs
        String? trimmomatic_version = read_QC_trim.trimmomatic_version
        String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
       
        # Read QC - fastp outputs
        String? fastp_version = read_QC_trim.fastp_version
        File? fastp_html_report = read_QC_trim.fastp_html_report
        
        # Read QC - bbduk outputs
        File? read1_clean = read_QC_trim.read1_clean
        File? read2_clean = read_QC_trim.read2_clean
        String? bbduk_docker = read_QC_trim.bbduk_docker
        
        # Read QC - cg pipeline outputs
        Float? r1_mean_q_raw = cg_pipeline_raw.r1_mean_q
        Float? r2_mean_q_raw = cg_pipeline_raw.r2_mean_q
        Float? combined_mean_q_raw = cg_pipeline_raw.combined_mean_q
        Float? r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength
        Float? r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength
        Float? combined_mean_readlength_raw = cg_pipeline_raw.combined_mean_readlength
        Float? r1_mean_q_clean = cg_pipeline_clean.r1_mean_q
        Float? r2_mean_q_clean = cg_pipeline_clean.r2_mean_q
        Float? combined_mean_q_clean = cg_pipeline_clean.combined_mean_q
        Float? r1_mean_readlength_clean = cg_pipeline_clean.r1_mean_readlength
        Float? r2_mean_readlength_clean = cg_pipeline_clean.r2_mean_readlength
        Float? combined_mean_readlength_clean = cg_pipeline_clean.combined_mean_readlength
        
        # Read QC - midas outputs
        String? midas_docker = read_QC_trim.midas_docker
        File? midas_report = read_QC_trim.midas_report
        String? midas_primary_genus = read_QC_trim.midas_primary_genus
        String? midas_secondary_genus = read_QC_trim.midas_secondary_genus
        Float? midas_secondary_genus_abundance = read_QC_trim.midas_secondary_genus_abundance
        Float? midas_secondary_genus_coverage = read_QC_trim.midas_secondary_genus_coverage
        
        # # Read QC - kraken outputs
        # String? kraken2_version = read_QC_trim.kraken_version
        # String? kraken2_report = read_QC_trim.kraken_report
        # String? kraken2_database = read_QC_trim.kraken_database
        # String? kraken_docker = read_QC_trim.kraken_docker

        # Assembly denovo outputs 
        File? assembly_fasta = assembly.assembly_fasta
        File? contigs_gfa = assembly.contigs_gfa
        File? filtered_contigs_metrics = assembly.filtered_contigs_metrics
        String? assembler = assembly.assembler_used
        String? assembler_version = assembly.assembler_version
        File? unicycler_log = assembly.unicycler_log
        Array[File]? unicycler_intermediate_graphs = assembly.unicycler_intermediate_graphs
        
        # Assembly QC - quast outputs
        File? quast_report = assembly.quast_report
        String? quast_version = assembly.quast_version
        Int? assembly_length = assembly.assembly_length
        Int? number_contigs = assembly.number_contigs
        Int? n50_value = assembly.n50_value
        Float? quast_gc_percent = assembly.quast_gc_percent
        
        # Assembly QC - cg pipeline outputs
        File? cg_pipeline_report_raw = cg_pipeline_raw.cg_pipeline_report
        String? cg_pipeline_docker = cg_pipeline_raw.cg_pipeline_docker
        Float? est_coverage_raw = cg_pipeline_raw.est_coverage
        File? cg_pipeline_report_clean = cg_pipeline_clean.cg_pipeline_report
        Float? est_coverage_clean = cg_pipeline_clean.est_coverage
        
        # Assembly QC - busco outputs
        String? busco_version = assembly.busco_version
        String? busco_docker = assembly.busco_docker
        String? busco_database = assembly.busco_database
        String? busco_results = assembly.busco_results
        File? busco_report = assembly.busco_report

        # Taxon ID - gambit outputs
        File? gambit_report = gambit.gambit_report_file
        File? gambit_closest_genomes = gambit.gambit_closest_genomes_file
        String? gambit_predicted_taxon = gambit.gambit_predicted_taxon
        String? gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
        String? gambit_version = gambit.gambit_version
        String? gambit_db_version = gambit.gambit_db_version
        String? gambit_docker = gambit.gambit_docker
        
		# ani-mummer outputs
        Float? ani_highest_percent = ani.ani_highest_percent
        Float? ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned
        File? ani_output_tsv = ani.ani_output_tsv
        String? ani_top_species_match = ani.ani_top_species_match
        String? ani_mummer_version = ani.ani_mummer_version
        String? ani_mummer_docker = ani.ani_docker
        
		# kmerfinder outputs
        String? kmerfinder_docker = kmerfinder.kmerfinder_docker
        File? kmerfinder_results_tsv = kmerfinder.kmerfinder_results_tsv
        String? kmerfinder_top_hit = kmerfinder.kmerfinder_top_hit
        String? kmerfinder_query_coverage = kmerfinder.kmerfinder_query_coverage
        String? kmerfinder_template_coverage = kmerfinder.kmerfinder_template_coverage
        String? kmerfinder_database = kmerfinder.kmerfinder_database
        
		# NCBI-AMRFinderPlus Outputs
        File? amrfinderplus_all_report = amrfinderplus_task.amrfinderplus_all_report
        File? amrfinderplus_amr_report = amrfinderplus_task.amrfinderplus_amr_report
        File? amrfinderplus_stress_report = amrfinderplus_task.amrfinderplus_stress_report
        File? amrfinderplus_virulence_report = amrfinderplus_task.amrfinderplus_virulence_report
        String? amrfinderplus_amr_core_genes = amrfinderplus_task.amrfinderplus_amr_core_genes
        String? amrfinderplus_amr_plus_genes = amrfinderplus_task.amrfinderplus_amr_plus_genes
        String? amrfinderplus_stress_genes = amrfinderplus_task.amrfinderplus_stress_genes
        String? amrfinderplus_virulence_genes = amrfinderplus_task.amrfinderplus_virulence_genes
        String? amrfinderplus_amr_classes = amrfinderplus_task.amrfinderplus_amr_classes
        String? amrfinderplus_amr_subclasses = amrfinderplus_task.amrfinderplus_amr_subclasses
        String? amrfinderplus_version = amrfinderplus_task.amrfinderplus_version
        String? amrfinderplus_db_version = amrfinderplus_task.amrfinderplus_db_version
        
        # NCBI-AMRFinderPlus Outputs for BETA-LACTAM genes
        String? amrfinderplus_amr_betalactam_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_genes
        String? amrfinderplus_amr_betalactam_betalactam_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_betalactam_genes
        String? amrfinderplus_amr_betalactam_carbapenem_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_carbapenem_genes
        String? amrfinderplus_amr_betalactam_cephalosporin_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_cephalosporin_genes
        String? amrfinderplus_amr_betalactam_cephalothin_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_cephalothin_genes
        String? amrfinderplus_amr_betalactam_methicillin_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_methicillin_genes
        
        # Resfinder Outputs
        File? resfinder_pheno_table = resfinder_task.resfinder_pheno_table
        File? resfinder_pheno_table_species = resfinder_task.resfinder_pheno_table_species
        File? resfinder_seqs = resfinder_task.resfinder_hit_in_genome_seq
        File? resfinder_results = resfinder_task.resfinder_results_tab
        File? resfinder_pointfinder_pheno_table = resfinder_task.pointfinder_pheno_table
        File? resfinder_pointfinder_results = resfinder_task.pointfinder_results
        String? resfinder_predicted_pheno_resistance = resfinder_task.resfinder_predicted_pheno_resistance
        String? resfinder_predicted_xdr_shigella = resfinder_task.resfinder_predicted_xdr_shigella
        String? resfinder_predicted_resistance_Amp = resfinder_task.resfinder_predicted_resistance_Amp
        String? resfinder_predicted_resistance_Azm = resfinder_task.resfinder_predicted_resistance_Azm
        String? resfinder_predicted_resistance_Axo = resfinder_task.resfinder_predicted_resistance_Axo
        String? resfinder_predicted_resistance_Cip = resfinder_task.resfinder_predicted_resistance_Cip
        String? resfinder_predicted_resistance_Smx = resfinder_task.resfinder_predicted_resistance_Smx
        String? resfinder_predicted_resistance_Tmp = resfinder_task.resfinder_predicted_resistance_Tmp
        String? resfinder_db_version = resfinder_task.resfinder_db_version
        String? resfinder_docker = resfinder_task.resfinder_docker
        
        # MLST Typing
        File? ts_mlst_results = ts_mlst.ts_mlst_results
        String? ts_mlst_predicted_st = ts_mlst.ts_mlst_predicted_st
        String? ts_mlst_pubmlst_scheme = ts_mlst.ts_mlst_pubmlst_scheme
        String? ts_mlst_allelic_profile = ts_mlst.ts_mlst_allelic_profile
        String? ts_mlst_version = ts_mlst.ts_mlst_version
        File? ts_mlst_novel_alleles = ts_mlst.ts_mlst_novel_alleles
        String? ts_mlst_docker = ts_mlst.ts_mlst_docker

        # AMR_Search
		File? amr_search_results = pathotype.amr_search_results
		File? amr_search_csv = pathotype.amr_results_csv
		File? amr_search_results_pdf = pathotype.amr_results_pdf
		String? amr_search_docker = pathotype.amr_search_docker
		String? amr_search_version = pathotype.amr_search_version
		
		# Ecoli Typing
		File? serotypefinder_report = pathotype.serotypefinder_report
		String? serotypefinder_docker = pathotype.serotypefinder_docker
		String? serotypefinder_serotype = pathotype.serotypefinder_serotype
		File? ectyper_results = pathotype.ectyper_results
		String? ectyper_version = pathotype.ectyper_version
		String? ectyper_predicted_serotype = pathotype.ectyper_predicted_serotype
		String? shigatyper_predicted_serotype = pathotype.shigatyper_predicted_serotype
		String? shigatyper_ipaB_presence_absence = pathotype.shigatyper_ipaB_presence_absence
		String? shigatyper_notes = pathotype.shigatyper_notes
		File? shigatyper_hits_tsv = pathotype.shigatyper_hits_tsv
		File? shigatyper_summary_tsv = pathotype.shigatyper_summary_tsv
		String? shigatyper_version = pathotype.shigatyper_version
		String? shigatyper_docker = pathotype.shigatyper_docker
		File? shigeifinder_report = pathotype.shigeifinder_report
		String? shigeifinder_docker = pathotype.shigeifinder_docker
		String? shigeifinder_version = pathotype.shigeifinder_version
		String? shigeifinder_ipaH_presence_absence = pathotype.shigeifinder_ipaH_presence_absence
		String? shigeifinder_num_virulence_plasmid_genes = pathotype.shigeifinder_num_virulence_plasmid_genes
		String? shigeifinder_cluster = pathotype.shigeifinder_cluster
		String? shigeifinder_serotype = pathotype.shigeifinder_serotype
		String? shigeifinder_O_antigen = pathotype.shigeifinder_O_antigen
		String? shigeifinder_H_antigen = pathotype.shigeifinder_H_antigen
		String? shigeifinder_notes = pathotype.shigeifinder_notes
		
		# ShigeiFinder outputs but for task that uses reads instead of assembly as input
		File? shigeifinder_report_reads = pathotype.shigeifinder_report
		String? shigeifinder_docker_reads = pathotype.shigeifinder_docker
		String? shigeifinder_version_reads = pathotype.shigeifinder_version
		String? shigeifinder_ipaH_presence_absence_reads = pathotype.shigeifinder_ipaH_presence_absence
		String? shigeifinder_num_virulence_plasmid_genes_reads = pathotype.shigeifinder_num_virulence_plasmid_genes
		String? shigeifinder_cluster_reads = pathotype.shigeifinder_cluster
		String? shigeifinder_serotype_reads = pathotype.shigeifinder_serotype
		String? shigeifinder_O_antigen_reads = pathotype.shigeifinder_O_antigen
		String? shigeifinder_H_antigen_reads = pathotype.shigeifinder_H_antigen
		String? shigeifinder_notes_reads = pathotype.shigeifinder_notes
		
		# E coli only typing
		File? virulencefinder_report_tsv = pathotype.virulencefinder_report_tsv
		String? virulencefinder_docker = pathotype.virulencefinder_docker
		String? virulencefinder_hits = pathotype.virulencefinder_hits
		
		# stxtyper 
		File? stxtyper_report = pathotype.stxtyper_report
		String? stxtyper_docker = pathotype.stxtyper_docker
		String? stxtyper_version = pathotype.stxtyper_version
		Int? stxtyper_num_hits = pathotype.stxtyper_num_hits
		String? stxtyper_all_hits = pathotype.stxtyper_all_hits
		String? stxtyper_complete_operons = pathotype.stxtyper_complete_operon_hits
		String? stxtyper_partial_hits = pathotype.stxtyper_partial_hits
		String? stxtyper_stx_frameshifts_or_internal_stop_hits =  pathotype.stxtyper_stx_frameshifts_or_internal_stop_hits
		String? stxtyper_novel_hits = pathotype.stxtyper_novel_hits
		String? stxtyper_extended_operons = pathotype.stxtyper_extended_operons
		String? stxtyper_ambiguous_hits = pathotype.stxtyper_ambiguous_hits
	
		# Shigella sonnei Typing
		File? sonneityping_mykrobe_report_csv = pathotype.sonneityping_mykrobe_report_csv
		File? sonneityping_mykrobe_report_json = pathotype.sonneityping_mykrobe_report_json
		File? sonneityping_final_report_tsv = pathotype.sonneityping_final_report_tsv
		String? sonneityping_mykrobe_version = pathotype.sonneityping_mykrobe_version
		String? sonneityping_mykrobe_docker = pathotype.sonneityping_mykrobe_docker
		String? sonneityping_species = pathotype.sonneityping_species
		String? sonneityping_final_genotype = pathotype.sonneityping_final_genotype
		String? sonneityping_genotype_confidence = pathotype.sonneityping_genotype_confidence
		String? sonneityping_genotype_name = pathotype.sonneityping_genotype_name
		
		# Listeria Typing
		File? lissero_results = pathotype.lissero_results
		String? lissero_version = pathotype.lissero_version
		String? lissero_serotype = pathotype.lissero_serotype
		
		# Pseudomonas Aeruginosa Typing
		String? pasty_serogroup = pathotype.pasty_serogroup
		Float? pasty_serogroup_coverage = pathotype.pasty_serogroup_coverage
		Int? pasty_serogroup_fragments = pathotype.pasty_serogroup_fragments
		File? pasty_summary_tsv = pathotype.pasty_summary_tsv
		File? pasty_blast_hits = pathotype.pasty_blast_hits
		File? pasty_all_serogroups = pathotype.pasty_all_serogroups
		String? pasty_version = pathotype.pasty_version
		String? pasty_docker = pathotype.pasty_docker
		String? pasty_comment = pathotype.pasty_comment
		
		# Salmonella Typing
		File? sistr_results = pathotype.sistr_results
		File? sistr_allele_json = pathotype.sistr_allele_json
		File? sistr_allele_fasta = pathotype.sistr_allele_fasta
		File? sistr_cgmlst = pathotype.sistr_cgmlst
		String? sistr_version = pathotype.sistr_version
		String? sistr_antigenic_formula = pathotype.sistr_antigenic_formula
		String? sistr_predicted_serotype = pathotype.sistr_predicted_serotype
		String? sistr_serogroup = pathotype.sistr_serogroup
		String? sistr_h1_antigens = pathotype.sistr_h1_antigens
		String? sistr_h2_antigens = pathotype.sistr_h2_antigens
		String? sistr_o_antigens = pathotype.sistr_o_antigens
		String? sistr_serotype_cgmlst = pathotype.sistr_serotype_cgmlst
		String? seqsero2_report = pathotype.seqsero2_report
		String? seqsero2_version = pathotype.seqsero2_version
		String? seqsero2_predicted_antigenic_profile = pathotype.seqsero2_predicted_antigenic_profile
		String? seqsero2_predicted_serotype = pathotype.seqsero2_predicted_serotype
		String? seqsero2_predicted_contamination = pathotype.seqsero2_predicted_contamination
		String? seqsero2_note = pathotype.seqsero2_note
		
		# Salmonella serotype Typhi Typing
		File? genotyphi_report_tsv = pathotype.genotyphi_report_tsv 
		File? genotyphi_mykrobe_json = pathotype.genotyphi_mykrobe_json
		String? genotyphi_version = pathotype.genotyphi_version
		String? genotyphi_species = pathotype.genotyphi_species
		Float? genotyphi_st_probes_percent_coverage = pathotype.genotyphi_st_probes_percent_coverage
		String? genotyphi_final_genotype = pathotype.genotyphi_final_genotype
		String? genotyphi_genotype_confidence = pathotype.genotyphi_genotype_confidence
		
		# Klebsiella Typing
		File? kleborate_output_file = pathotype.kleborate_output_file
		String? kleborate_version = pathotype.kleborate_version
		String? kleborate_docker = pathotype.kleborate_docker
		String? kleborate_key_resistance_genes = pathotype.kleborate_key_resistance_genes
		String? kleborate_genomic_resistance_mutations = pathotype.kleborate_genomic_resistance_mutations
		String? kleborate_mlst_sequence_type = pathotype.kleborate_mlst_sequence_type
		String? kleborate_klocus = pathotype.kleborate_klocus
		String? kleborate_ktype = pathotype.kleborate_ktype
		String? kleborate_olocus = pathotype.kleborate_olocus
		String? kleborate_otype = pathotype.kleborate_otype
		String? kleborate_klocus_confidence = pathotype.kleborate_klocus_confidence
		String? kleborate_olocus_confidence = pathotype.kleborate_olocus_confidence
		String? kleborate_virulence_score = pathotype.kleborate_virulence_score
		String? kleborate_resistance_score = pathotype.kleborate_resistance_score
		
		# Neisseria gonorrhoeae Typing
		File? ngmaster_tsv = pathotype.ngmaster_tsv
		String? ngmaster_version = pathotype.ngmaster_version
		String? ngmaster_ngmast_sequence_type = pathotype.ngmaster_ngmast_sequence_type
		String? ngmaster_ngmast_porB_allele = pathotype.ngmaster_ngmast_porB_allele
		String? ngmaster_ngmast_tbpB_allele = pathotype.ngmaster_ngmast_tbpB_allele
		String? ngmaster_ngstar_sequence_type = pathotype.ngmaster_ngstar_sequence_type
		String? ngmaster_ngstar_penA_allele = pathotype.ngmaster_ngstar_penA_allele
		String? ngmaster_ngstar_mtrR_allele = pathotype.ngmaster_ngstar_mtrR_allele
		String? ngmaster_ngstar_porB_allele = pathotype.ngmaster_ngstar_porB_allele
		String? ngmaster_ngstar_ponA_allele = pathotype.ngmaster_ngstar_ponA_allele
		String? ngmaster_ngstar_gyrA_allele = pathotype.ngmaster_ngstar_gyrA_allele
		String? ngmaster_ngstar_parC_allele = pathotype.ngmaster_ngstar_parC_allele
		String? ngmaster_ngstar_23S_allele = pathotype.ngmaster_ngstar_23S_allele
		
		# Neisseria meningitidis Typing
		File? meningotype_tsv = pathotype.meningotype_tsv
		String? meningotype_version = pathotype.meningotype_version
		String? meningotype_serogroup = pathotype.meningotype_serogroup
		String? meningotype_PorA = pathotype.meningotype_PorA
		String? meningotype_FetA = pathotype.meningotype_FetA
		String? meningotype_PorB = pathotype.meningotype_PorB
		String? meningotype_fHbp = pathotype.meningotype_fHbp
		String? meningotype_NHBA = pathotype.meningotype_NHBA
		String? meningotype_NadA = pathotype.meningotype_NadA
		String? meningotype_BAST = pathotype.meningotype_BAST
		
		# Acinetobacter Typing
		File? kaptive_output_file_k = pathotype.kaptive_output_file_k
		File? kaptive_output_file_oc = pathotype.kaptive_output_file_oc
		String? kaptive_version = pathotype.kaptive_version
		String? kaptive_k_locus = pathotype.kaptive_k_match
		String? kaptive_k_type = pathotype.kaptive_k_type
		String? kaptive_kl_confidence = pathotype.kaptive_k_confidence
		String? kaptive_oc_locus = pathotype.kaptive_oc_match
		String? kaptive_ocl_confidence = pathotype.kaptive_oc_confidence
		File? abricate_abaum_plasmid_tsv = pathotype.abricate_abaum_results
		String? abricate_abaum_plasmid_type_genes = pathotype.abricate_abaum_genes
		String? abricate_abaum_database = pathotype.abricate_abaum_database
		String? abricate_abaum_version = pathotype.abricate_abaum_version
		String? abricate_abaum_docker = pathotype.abricate_abaum_docker
		
		# Mycobacterium Typing
		File? tbprofiler_output_file = pathotype.tbprofiler_output_file
		File? tbprofiler_output_bam = pathotype.tbprofiler_output_bam
		File? tbprofiler_output_bai = pathotype.tbprofiler_output_bai
		File? tbprofiler_output_vcf = pathotype.tbprofiler_output_vcf
		String? tbprofiler_version = pathotype.tbprofiler_version
		String? tbprofiler_main_lineage = pathotype.tbprofiler_main_lineage
		String? tbprofiler_sub_lineage = pathotype.tbprofiler_sub_lineage
		String? tbprofiler_dr_type = pathotype.tbprofiler_dr_type
		String? tbprofiler_resistance_genes = pathotype.tbprofiler_resistance_genes
		Float? tbprofiler_median_depth = pathotype.tbprofiler_median_depth
		Float? tbprofiler_pct_reads_mapped = pathotype.tbprofiler_pct_reads_mapped
		String? tbp_parser_version = pathotype.tbp_parser_version
		String? tbp_parser_docker = pathotype.tbp_parser_docker
		File? tbp_parser_lims_report_csv = pathotype.tbp_parser_lims_report_csv
		File? tbp_parser_looker_report_csv = pathotype.tbp_parser_looker_report_csv
		File? tbp_parser_laboratorian_report_csv = pathotype.tbp_parser_laboratorian_report_csv
		File? tbp_parser_coverage_report = pathotype.tbp_parser_coverage_report
		Float? tbp_parser_genome_percent_coverage = pathotype.tbp_parser_genome_percent_coverage
		Float? tbp_parser_average_genome_depth = pathotype.tbp_parser_average_genome_depth
		File? clockwork_decontaminated_read1 = pathotype.clockwork_cleaned_read1
		File? clockwork_decontaminated_read2 = pathotype.clockwork_cleaned_read2
		
		# Legionella pneumophila typing
		File? legsta_results = pathotype.legsta_results
		String? legsta_predicted_sbt = pathotype.legsta_predicted_sbt
		String? legsta_version = pathotype.legsta_version
		
		# Staphylococcus aureus
		File? spatyper_tsv = pathotype.spatyper_tsv
		String? spatyper_docker = pathotype.spatyper_docker
		String? spatyper_repeats = pathotype.spatyper_repeats
		String? spatyper_type = pathotype.spatyper_type
		String? spatyper_version = pathotype.spatyper_version
		File? staphopiasccmec_results_tsv = pathotype.staphopiasccmec_results_tsv
		File? staphopiasccmec_hamming_distance_tsv = pathotype.staphopiasccmec_hamming_distance_tsv
		String? staphopiasccmec_types_and_mecA_presence = pathotype.staphopiasccmec_types_and_mecA_presence
		String? staphopiasccmec_version = pathotype.staphopiasccmec_version
		String? staphopiasccmec_docker = pathotype.staphopiasccmec_docker
		File? agrvate_summary = pathotype.agrvate_summary
		File? agrvate_results = pathotype.agrvate_results
		String? agrvate_agr_group = pathotype.agrvate_agr_group
		String? agrvate_agr_match_score = pathotype.agrvate_agr_match_score
		String? agrvate_agr_canonical = pathotype.agrvate_agr_canonical
		String? agrvate_agr_multiple = pathotype.agrvate_agr_multiple
		String? agrvate_agr_num_frameshifts = pathotype.agrvate_agr_num_frameshifts
		String? agrvate_version = pathotype.agrvate_version
		String? agrvate_docker = pathotype.agrvate_docker
		
		# Streptococcus pneumoniae Typing
		String? pbptyper_predicted_1A_2B_2X = pathotype.pbptyper_predicted_1A_2B_2X
		File? pbptyper_pbptype_predicted_tsv = pathotype.pbptyper_pbptype_predicted_tsv
		String? pbptyper_version = pathotype.pbptyper_version
		String? pbptyper_docker = pathotype.pbptyper_docker
		String? poppunk_gps_cluster = pathotype.poppunk_gps_cluster
		File? poppunk_gps_external_cluster_csv = pathotype.poppunk_gps_external_cluster_csv
		String? poppunk_GPS_db_version = pathotype.poppunk_GPS_db_version
		String? poppunk_version = pathotype.poppunk_version
		String? poppunk_docker = pathotype.poppunk_docker
		String? seroba_version = pathotype.seroba_version
		String? seroba_docker = pathotype.seroba_docker
		String? seroba_serotype = pathotype.seroba_serotype
		String? seroba_ariba_serotype = pathotype.seroba_ariba_serotype
		String? seroba_ariba_identity = pathotype.seroba_ariba_identity
		File? seroba_details = pathotype.seroba_details
		
		# Streptococcus pyogenes Typing
		String? emmtyper_emm_type = pathotype.emmtyper_emm_type
		File? emmtyper_results_tsv = pathotype.emmtyper_results_tsv
		String? emmtyper_version = pathotype.emmtyper_version
		String? emmtyper_docker = pathotype.emmtyper_docker
		String? emmtypingtool_emm_type = pathotype.emmtypingtool_emm_type
		File? emmtypingtool_results_xml = pathotype.emmtypingtool_results_xml
		String? emmtypingtool_version = pathotype.emmtypingtool_version
		String? emmtypingtool_docker = pathotype.emmtypingtool_docker
		
		# Haemophilus influenzae Typing
		String? hicap_serotype = pathotype.hicap_serotype
		String? hicap_genes = pathotype.hicap_genes
		File? hicap_results_tsv = pathotype.hicap_results_tsv
		String? hicap_version = pathotype.hicap_version
		String? hicap_docker = pathotype.hicap_docker
	
		# Vibrio Typing
		File? srst2_vibrio_detailed_tsv = pathotype.srst2_vibrio_detailed_tsv
		String? srst2_vibrio_version = pathotype.srst2_vibrio_version
		String? srst2_vibrio_docker = pathotype.srst2_vibrio_docker
		String? srst2_vibrio_database = pathotype.srst2_vibrio_database
		String? srst2_vibrio_ctxA = pathotype.srst2_vibrio_ctxA
		String? srst2_vibrio_ompW = pathotype.srst2_vibrio_ompW
		String? srst2_vibrio_toxR = pathotype.srst2_vibrio_toxR
		String? srst2_vibrio_biotype = pathotype.srst2_vibrio_biotype
		String? srst2_vibrio_serogroup = pathotype.srst2_vibrio_serogroup
		File? abricate_vibrio_detailed_tsv = pathotype.abricate_vibrio_detailed_tsv
		String? abricate_vibrio_database = pathotype.abricate_vibrio_database
		String? abricate_vibrio_docker = pathotype.abricate_vibrio_docker
		String? abricate_vibrio_version = pathotype.abricate_vibrio_version
		String? abricate_vibrio_ctxA = pathotype.abricate_vibrio_ctxA
		String? abricate_vibrio_ompW = pathotype.abricate_vibrio_ompW
		String? abricate_vibrio_toxR = pathotype.abricate_vibrio_toxR
		String? abricate_vibrio_biotype = pathotype.abricate_vibrio_biotype
		String? abricate_vibrio_serogroup = pathotype.abricate_vibrio_serogroup
		File? vibecheck_lineage_report = pathotype.vibecheck_lineage_report
		String? vibecheck_top_lineage = pathotype.vibecheck_top_lineage
		Float? vibecheck_confidence = pathotype.vibecheck_confidence
		String? vibecheck_classification_notes = pathotype.vibecheck_classification_notes
		String? vibecheck_version = pathotype.vibecheck_version
		String? vibecheck_docker = pathotype.vibecheck_docker

    }

}
