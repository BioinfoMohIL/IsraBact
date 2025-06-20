version 1.0

import "../../tasks/assembly/task_unicycler.wdl" as task_unicycler
import "../../tasks/alignment/task_bwa.wdl" as task_bwa
import "../../tasks/qc/read_filtering/task_pilon.wdl" as task_pilon
import "../../tasks/qc/advanced_metrics/task_busco.wdl" as busco_task
import "../../tasks/qc/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/qc/read_filtering/task_filter_contigs.wdl" as task_filter_contigs
import "./wf_read_qc_trim_pe.wdl" as read_QC


workflow assembly_long_reads {
	meta {
		description: "De-novo genome assembly for long reads"
		author: "David Maimoun"
	}

	input {
		File read1
		File read2
		String samplename
		String assembler = "unicycler" # Options: unicycler
		String seq_method = 'Illumina'
	
		# trimming parameters
        Int trim_min_length = 75
        Int trim_quality_min_score = 20
        Int trim_window_size = 10

		# Optional parameters for unicycler
		File? illumina_unpaired_fq
		File? long_reads

		Boolean call_pilon = false
		String? assembler_options # Extra assembler options
		Boolean run_filter_contigs = false # Default: Filter contigs after assembly

		# Optional parameters for bwa
		Int? bwa_cpu
		Int? bwa_memory
		Int? bwa_disk_size
		String? bwa_docker
		
		# Optional parameters for pilon
		Int? pilon_cpu
		Int? pilon_memory
		Int? pilon_disk_size
		String? pilon_docker
		Int? pilon_min_mapping_quality = 60 # Shovill default
		Int? pilon_min_base_quality = 3 # Shovill default
		Float? pilon_min_depth = 0.25 # Shovill default
		String? pilon_fix = "bases" # Options: all, snps, indels, gaps
		
		# Optional parameters for filtering 
		Int filter_contigs_min_length = 200 # Default we set before
		Float filter_contigs_min_coverage = 2.0 # Default we set before
		Boolean filter_contigs_skip_length_filter = false
		Boolean filter_contigs_skip_coverage_filter = false
		Boolean filter_contigs_skip_homopolymer_filter = false
		Int? filter_contigs_cpu
		Int? filter_contigs_memory
		Int? filter_contigs_disk_size
		String? filter_contigs_docker
		
	}

	call read_QC.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1 = select_first([read1]),
        read2 = read2,
        trim_min_length = trim_min_length,
        trim_quality_min_score = trim_quality_min_score,
        trim_window_size = trim_window_size,
        workflow_series = "israbact_pe"
    }
	
	call task_unicycler.unicycler {
		input:
			read1 = read1,
			read2 = read2,
			samplename = samplename,
			illumina_unpaired = illumina_unpaired_fq,
			long_reads = long_reads
	}

	if (call_pilon) {
		call task_bwa.bwa {
		input:
			read1 = read1,
			read2 = read2,
			samplename = samplename,
			reference_genome = select_first([unicycler.assembly_fasta]),
			cpu = bwa_cpu,
			memory = bwa_memory,
			disk_size = bwa_disk_size,
			docker = bwa_docker
		}

		call task_pilon.pilon {
			input:
				assembly = select_first([unicycler.assembly_fasta]),
				bam = bwa.sorted_bam,
				bai = bwa.sorted_bai,
				samplename = samplename,
				min_mapping_quality = pilon_min_mapping_quality,
				min_base_quality = pilon_min_base_quality,
				min_depth = pilon_min_depth,
				fix = pilon_fix,
				cpu = pilon_cpu,
				memory = pilon_memory,
				disk_size = pilon_disk_size,
				docker = pilon_docker
		}
	}

	if (run_filter_contigs) {
		call task_filter_contigs.filter_contigs {
			input:
				samplename = samplename,
				assembly_fasta = select_first([unicycler.assembly_fasta]),
				min_length = filter_contigs_min_length,
				min_coverage = filter_contigs_min_coverage,
				skip_length_filter = filter_contigs_skip_length_filter,
				skip_coverage_filter = filter_contigs_skip_coverage_filter,
				skip_homopolymer_filter = filter_contigs_skip_homopolymer_filter,
				cpu = filter_contigs_cpu,
				memory = filter_contigs_memory,
				disk_size = filter_contigs_disk_size,
				docker = filter_contigs_docker
		}
	}

	File assembly_file = select_first([unicycler.assembly_fasta])

	call quast_task.quast {
		input:
			assembly = assembly_file,
			samplename = samplename
	}

	call busco_task.busco {
		input:
			assembly = assembly_file,
			samplename = samplename
	}

	
	output {
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
        
    
		# Assembler outputs
		File lr_assembly_fasta = assembly_file
		String lr_assembler_used = assembler
		File? lr_contigs_gfa = unicycler.assembly_gfa
		File? lr_unicycler_log = unicycler.unicycler_log
		Array[File]? lr_unicycler_intermediate_graphs = unicycler.unicycler_intermediate_graphs
		String? lr_assembler_version = unicycler.unicycler_version

		File? lr_filtered_contigs_metrics = filter_contigs.assembly_filtering_metrics
		File? lr_pilon_changes = pilon.changes
		File? lr_pilon_vcf = pilon.vcf

		# Assembly QC - quast outputs
		File? lr_quast_report = quast.quast_report
		String? lr_quast_version = quast.version
		Int? lr_assembly_length = quast.genome_length
		Int? lr_number_contigs = quast.number_contigs
		Int? lr_n50_value = quast.n50_value
		Float? lr_quast_gc_percent = quast.gc_percent

		# Assembly QC - busco outputs
		String? lr_busco_version = busco.busco_version
		String? lr_busco_docker = busco.busco_docker
		String? lr_busco_database = busco.busco_database
		String? lr_busco_results = busco.busco_results
		File? lr_busco_report = busco.busco_report
	}

}
   
