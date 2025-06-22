version 1.0

import "../../../tasks/species_typing/neisseria/neisseria_typing.wdl" as nm_typing
import "../../../tasks/qc/read/task_cg_pipeline.wdl" as cg
import "../../utilities/wf_assembly.wdl" as deno_assembly
import "../../utilities/wf_read_qc_trim_pe.wdl" as read_QC
import "../../utilities/file_handling/wf_concatenate_illumina_lanes.wdl" as concatenate_lanes_workflow

workflow neisseria_typing {
	meta {
		description: "N. Meningitidis De-novo genome assembly, QC, and typing of paired-end NGS data"
		author: "David Maimoun (The Codon Bleu)"
		email: "thecodonbleu@outlook.com"
	}

	input {
		String samplename
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

		# trimming parameters
        Int trim_min_length = 75
        Int trim_quality_min_score = 20
        Int trim_window_size = 10

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
			assembler = 'spades'
	}

	# call cg.cg_pipeline as cg_pipeline_raw {
	# 	input:
	# 		read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
	# 		read2 = select_first([concatenate_illumina_lanes.read2_concatenated, read2]),
	# 		samplename = samplename,
	# 		genome_length = select_first([genome_length, assembly.assembly_length])
	# }

	# call cg.cg_pipeline as cg_pipeline_clean {
	# 	input:
	# 		read1 = read_QC_trim.read1_clean,
	# 		read2 = read_QC_trim.read2_clean,
	# 		samplename = samplename,
	# 		genome_length = select_first([genome_length, assembly.assembly_length])
	# }


    call nm_typing.neisseria_typing as typing {
		input:
			assembly = assembly.assembly_fasta
    }

    call nm_typing.serogrouping as serogrouping {
        input: 
          assembly = assembly.assembly_fasta,
          samplename = samplename,
    }

    output {
        # QC - outputs
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
        
  
       	# Assembly denovo outputs 
		String assembler_used = assembly.assembler_used
		String seq_platform = assembly.seq_platform

        File? assembly_fasta = assembly.assembly_fasta
        File? contigs_gfa = assembly.contigs_gfa
        File? filtered_contigs_metrics = assembly.filtered_contigs_metrics
        String? assembler = assembly.assembler_used
        String? assembler_version = assembly.assembler_version
        
        # Assembly QC - quast outputs
        File? quast_report = assembly.quast_report
        String? quast_version = assembly.quast_version
        Int? assembly_length = assembly.assembly_length
        Int? number_contigs = assembly.number_contigs
        Int? n50_value = assembly.n50_value
        Float? quast_gc_percent = assembly.quast_gc_percent
           
        # Assembly QC - busco outputs
        String? busco_version = assembly.busco_version
        File? busco_report = assembly.busco_report
        String? busco_docker = assembly.busco_docker
        String? busco_database = assembly.busco_database
        String? busco_results = assembly.busco_results

         # Assembly QC - cg pipeline outputs
        # File? cg_pipeline_report_raw = cg_pipeline_raw.cg_pipeline_report
        # String? cg_pipeline_docker = cg_pipeline_raw.cg_pipeline_docker
        # Float? est_coverage_raw = cg_pipeline_raw.est_coverage
        # File? cg_pipeline_report_clean = cg_pipeline_clean.cg_pipeline_report
        # Float? est_coverage_clean = cg_pipeline_clean.est_coverage
    
        # Typing
        # File typing_results = typing.typing_results
        String neisseria_typing_version = typing.version
        String fHbp_peptide = typing.fHbp_peptide
        String NHBA_peptide = typing.NHBA_peptide
        String NadA_peptide = typing.NadA_peptide
        String PorA_VR1 = typing.PorA_VR1 
        String PorA_VR2 = typing.PorA_VR2 
        String abcZ = typing.abcZ
        String adk = typing.adk
        String aroE = typing.aroE
        String gdh = typing.gdh
        String fumC = typing.fumC
        String pdhC = typing.pdhC
        String pgm = typing.pgm
        String st = typing.st 
        String FetA_VR = typing.FetA_VR
        String bast_type = typing.bast_type
        String clonal_complex = typing.clonal_complex
        String mendevar_bexsero_reactivity  = typing.bexsero_cross_reactivity      
        String mendevar_trumenba_reactivity = typing.trumenba_cross_reactivity

        String serogroup = serogrouping.serogroup
        String genogroup = serogrouping.genogroup
  
    }
   
    
}


