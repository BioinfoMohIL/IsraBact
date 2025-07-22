version 1.0

import "../../tasks/assembly/task_fq2dna.wdl" as task_fq2dna

workflow wf_fq2dna {
  input {
    File read1
    File read2
    String species
    String organism = "B"     # "B(Bacteria), P(Prokariote),  E(Eukaryote), S(Standard), V(Virus)"
    String alien_tag = "AUTO"
  }

  meta {
		description: "De-novo genome assembly, QC your NGS reads from fq2DNA (Alexis Criscuolo - Pasteur Institut) "
		author: "David Maimoun"
    organization: "MOH"
	

  call task_fq2dna.fq2dna_run {
    input:
        read1 = read1,
        read2 = read2,
        species = species,
        organism = organism,
        alien_tag = alien_tag
    }

  
  output {
    String fq2dna_version           = fq2dna_run.version
    File fq2dna_assembly_fasta      = fq2dna_run.assembly_fasta
    File fq2dna_metrics_zip         = fq2dna_run.metrics_zip
    File fq2dna_selected_scaffolds  = fq2dna_run.selected_scaffolds
    File fq2dna_selected_contigs    = fq2dna_run.selected_contigs
    File fq2dna_scaffolding_info    = fq2dna_run.scaffolding_info
  }
}

