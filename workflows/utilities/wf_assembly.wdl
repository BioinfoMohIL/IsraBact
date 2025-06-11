version 1.0

import "../../tasks/qc/task_trimming.wdl" as trimming
import "../../tasks/qc/task_assembly_qc.wdl" as assembly_qc
import "../../tasks/qc/task_reads_qc.wdl" as reads_qc
import "../../tasks/qc/task_cg_pipeline.wdl" as cg
import "../../tasks/assembly/task_assembly.wdl" as assemble_reads


workflow assembly {
  meta {
    description: "De-novo genome assembly, QC your NGS reads"
    author: "David Maimoun"
    organization: "MOH Jerusalem"
  }

  input {
    File read1
    File read2
    String samplename
    String seq_method = "ILLUMINA"

    Int? genome_length
 
    Int trim_min_length = 75
    Int trim_quality_min_score = 20
    Int trim_window_size = 10
    
  }
  
    call trimming.trimmomatic as trimmomatic {
      input:
        samplename = samplename,
        read1 = read1,
        read2 = read2 ,
        trimmomatic_window_size = trim_window_size,
        trimmomatic_quality_trim_score = trim_quality_min_score,
        trimmomatic_min_length = trim_min_length,
    }
  
    call reads_qc.fastqc as fastqc {
      input:
        read1 = trimmomatic.read1_trimmed,
        read2 = trimmomatic.read2_trimmed
    }
    
    call assemble_reads.shovill as shovill {
        input:
            samplename = samplename,
            read1 = trimmomatic.read1_trimmed,
            read2 = trimmomatic.read2_trimmed
      }

    call assembly_qc.quast as quast {
        input:
          assembly = shovill.assembly_fasta,
          samplename = samplename
    }

    call cg.cg_pipeline as cg_pipeline_raw {
        input:
          read1 = trimmomatic.read1_trimmed,
          read2 = trimmomatic.read2_trimmed,
          samplename = samplename,
          genome_length = select_first([genome_length, quast.genome_length])
      
    }

    call cg.cg_pipeline as cg_pipeline_clean {
      input:
        read1 = trimmomatic.read1_trimmed,
        read2 = trimmomatic.read2_trimmed,
        samplename = samplename,
        genome_length = select_first([genome_length, quast.genome_length])
    }
    
    call reads_qc.busco as busco {
      input:
        assembly = shovill.assembly_fasta,
        samplename = samplename
    }


    output {
        String seq_platform = seq_method
    
        # Trimmomatic outputs
        String trimmomatic_version = trimmomatic.version
        String trimmomatic_docker = trimmomatic.trimmomatic_docker
        
        # Read QC - FastQC outputs
        File read1_fastqc_html = fastqc.read1_fastqc_html
        File read2_fastqc_html = fastqc.read2_fastqc_html
        String version = fastqc.version
        String fastqc_docker = fastqc.fastqc_docker
        File read1_fastqc_zip  = fastqc.read1_fastqc_zip
        File read2_fastqc_zip  = fastqc.read2_fastqc_zip
        Int read1_seq = fastqc.read1_seq
        Int read2_seq = fastqc.read2_seq

        # Read QC - cg pipeline outputs
        Float r1_mean_q_raw = cg_pipeline_raw.r1_mean_q
        Float r2_mean_q_raw = cg_pipeline_raw.r2_mean_q
        Float combined_mean_q_raw = cg_pipeline_raw.combined_mean_q
        Float r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength
        Float r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength
        Float combined_mean_readlength_raw = cg_pipeline_raw.combined_mean_readlength
        Float r1_mean_q_clean = cg_pipeline_clean.r1_mean_q
        Float r2_mean_q_clean = cg_pipeline_clean.r2_mean_q
        Float combined_mean_q_clean = cg_pipeline_clean.combined_mean_q
        Float r1_mean_readlength_clean = cg_pipeline_clean.r1_mean_readlength
        Float r2_mean_readlength_clean = cg_pipeline_clean.r2_mean_readlength
        Float combined_mean_readlength_clean = cg_pipeline_clean.combined_mean_readlength
  
        # Assembly - shovill outputs 
        File assembly_fasta = shovill.assembly_fasta
        String shovill_version = shovill.shovill_version
        File? contigs_gfa = shovill.contigs_gfa
        File? contigs_fastg = shovill.contigs_fastg
        File? contigs_lastgraph = shovill.contigs_lastgraph
        
        # Assembly QC - quast outputs
        File quast_report = quast.quast_report
        String quast_version = quast.version
        Int assembly_length = quast.genome_length
        Int number_contigs = quast.number_contigs
        Int n50_value = quast.n50_value
        Float quast_gc_percent = quast.gc_percent
           
        # Assembly QC - busco outputs
        String busco_version = busco.busco_version
        File busco_report = busco.busco_report
        String busco_docker = busco.busco_docker
        String busco_database = busco.busco_database
        String busco_results = busco.busco_results

        # Assembly QC - cg pipeline outputs
        File cg_pipeline_report_raw = cg_pipeline_raw.cg_pipeline_report
        Float est_coverage_raw = cg_pipeline_raw.est_coverage
        File cg_pipeline_report_clean = cg_pipeline_clean.cg_pipeline_report
        Float est_coverage_clean = cg_pipeline_clean.est_coverage
        String cg_pipeline_docker = cg_pipeline_raw.cg_pipeline_docker
    }
   
    
}

