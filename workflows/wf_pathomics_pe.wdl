version 1.0

import "../tasks/task_versioning.wdl" as versioning
import "../tasks/qc/validation/validate_reads.wdl" as validate_reads
import "../tasks/qc/read/task_cg_pipeline.wdl" as cg
import "../tasks/qc/assembly/task_busco.wdl" as busco_task
import "../tasks/qc/assembly/task_quast.wdl" as quast_task
import "./utilities/wf_read_qc_trim_pe.wdl" as read_QC
import "./utilities/wf_assembly.wdl" as deno_assembly
import "./utilities/file_handling/wf_concatenate_illumina_lanes.wdl" as concatenate_lanes_workflow


workflow wf_pathomics_pe {
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
        workflow_series = "pathomics"
    }

    call deno_assembly.assembly {
        input:
            samplename = samplename,
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            illumina_unpaired_fq = illumina_unpaired_fq,
            long_reads = long_reads
    }

    call quast_task.quast {
        input:
          assembly = assembly.assembly_fasta,
          samplename = samplename
      }

      call cg.cg_pipeline as cg_pipeline_raw {
        input:
          read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
          read2 = select_first([concatenate_illumina_lanes.read2_concatenated, read2]),
          samplename = samplename,
          genome_length = select_first([genome_length, quast.genome_length])
      }

      call cg.cg_pipeline as cg_pipeline_clean {
        input:
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          samplename = samplename,
          genome_length = select_first([genome_length, quast.genome_length])
      }

      call busco_task.busco {
        input:
          assembly = assembly.assembly_fasta,
          samplename = samplename
      }

    output {
        String omics_du_chef_pe_version = version_capture.version
        String omics_du_chef_pe_date    = version_capture.date
        
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
        File? quast_report = quast.quast_report
        String? quast_version = quast.version
        Int? assembly_length = quast.genome_length
        Int? number_contigs = quast.number_contigs
        Int? n50_value = quast.n50_value
        Float? quast_gc_percent = quast.gc_percent
        
        # Assembly QC - cg pipeline outputs
        File? cg_pipeline_report_raw = cg_pipeline_raw.cg_pipeline_report
        String? cg_pipeline_docker = cg_pipeline_raw.cg_pipeline_docker
        Float? est_coverage_raw = cg_pipeline_raw.est_coverage
        File? cg_pipeline_report_clean = cg_pipeline_clean.cg_pipeline_report
        Float? est_coverage_clean = cg_pipeline_clean.est_coverage
        
        # Assembly QC - busco outputs
        String? busco_version = busco.busco_version
        String? busco_docker = busco.busco_docker
        String? busco_database = busco.busco_database
        String? busco_results = busco.busco_results
        File? busco_report = busco.busco_report
    }

}
