version 1.0

import "../tasks/typing/neisseria/neisseria_typing.wdl" as nm_typing
import "./utilities/wf_assembly.wdl" as assemble_reads

workflow neisseria_typing {
  meta {
    description: "N. Meningitidis De-novo genome assembly, QC, and typing of paired-end NGS data"
    author: "David Maimoun"
    organization: "MOH Jerusalem"
  }

  input {
    String samplename
    File read1
    File read2  
  }
  

  
  call assemble_reads.assembly {
    input:
        samplename = samplename,
        read1 = read1,
        read2 = read2
  }


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
        String seq_platform = assembly.seq_platform
    
        # Trimmomatic outputs
        String trimmomatic_version = assembly.trimmomatic_version
        String trimmomatic_docker = assembly.trimmomatic_docker 
        
        # Read QC - FastQC outputs
        File fastqc_read1_html = assembly.read1_fastqc_html 
        File fastqc_read2_html = assembly.read2_fastqc_html
        String fastqc_version = assembly.version
        String fastqc_docker = assembly.fastqc_docker
        File fastqc_read1_zip  = assembly.read1_fastqc_zip
        File fastqc_read2_zip  = assembly.read2_fastqc_zip
      

        # Read QC - cg pipeline outputs
        Float r1_mean_q_raw = assembly.r1_mean_q_raw
        Float r2_mean_q_raw = assembly.r2_mean_q_raw 
        Float combined_mean_q_raw = assembly.combined_mean_q_raw
        Float r1_mean_readlength_raw = assembly.r1_mean_readlength_raw
        Float r2_mean_readlength_raw = assembly.r2_mean_readlength_raw
        Float combined_mean_readlength_raw = assembly.combined_mean_readlength_raw
        Float r1_mean_q_clean = assembly.r1_mean_q_clean 
        Float r2_mean_q_clean = assembly.r2_mean_q_clean
        Float combined_mean_q_clean = assembly.combined_mean_q_clean 
        Float r1_mean_readlength_clean = assembly.r1_mean_readlength_clean
        Float r2_mean_readlength_clean = assembly.r2_mean_readlength_clean
        Float combined_mean_readlength_clean = assembly.combined_mean_readlength_clean
  
        # Assembly - shovill outputs 
        File assembly_fasta = assembly.assembly_fasta
        String shovill_version = assembly.shovill_version
        File? contigs_gfa = assembly.contigs_gfa
        File? contigs_fastg = assembly.contigs_fastg
        File? contigs_lastgraph = assembly.contigs_lastgraph
        
        # Assembly QC - quast outputs
        File quast_report = assembly.quast_report
        String quast_version = assembly.quast_version
        Int assembly_length = assembly.assembly_length
        Int number_contigs = assembly.number_contigs
        Int n50_value = assembly.n50_value
        Float quast_gc_percent = assembly.quast_gc_percent
           
        # Assembly QC - busco outputs
        String busco_version = assembly.busco_version
        File busco_report = assembly.busco_report
        String busco_docker = assembly.busco_docker
        String busco_database = assembly.busco_database
        String busco_results = assembly.busco_results

        # Assembly QC - cg pipeline outputs
        File cg_pipeline_report_raw = assembly.cg_pipeline_report_raw
        Float est_coverage_raw = assembly.est_coverage_raw
        File cg_pipeline_report_clean = assembly.cg_pipeline_report_clean
        Float est_coverage_clean = assembly.est_coverage_clean
        String cg_pipeline_docker = assembly.cg_pipeline_docker
    

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


