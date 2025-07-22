IsraBact WDL Workflows
=======================================

De-novo genome assembly, taxonomic ID, genes and species typing, alleles comparison, QC of paired-end bacterial NGS data, based on Theiagen - Public Health Bioinformatics workflows


This repository provides WDL (Workflow Description Language) scripts for various bioinformatics workflows related to bacterial genomics. The workflows are designed to be run using MiniWDL or Cromwell with Docker support.

Each workflow performs a specific task related to genome analysis, quality control, or virulence factor detection. Below are the detailed descriptions for each script.


-------------------------------------------------------------------------------
1. Workflow: metagenomics/wf_species_detection_bs_reads.wdl
-------------------------------------------------------------------------------

Description:
  Check whether the sequenced samples in Basespace Illumina correspond to the expected species using Kraken2.

Inputs:
  - **[String]** api_server: "Illumina Basepace API"
  - *[String]* access_token: "Illumina Basespace access token"
  - [String] basespace_collection_id: "Samplename in the Basespace platform (for ex, EC001, NM005 , = our entity id)"
  - [String] sample_prefix: "Optional, to fetch only specific species according to your samplename prefix, for ex 'EC' for ecoli (useful for testing)"
    
Outputs:
  - **[String]** version.
  - **[File]** Reads list file: "reads_list.txt".txt.
  - **[Array[String]]** list of samples name.
  - **[File]** Species detected table: species_detected.csv. Display the samplename (Sample), the species detected by Kraken (Detected), and + if detected matchs to expected (via the prefix, for ex, if samplename EC001, expected ecoli)


-------------------------------------------------------------------------------
2. Workflow: virulence/wf_virulence_finder.wdl
-------------------------------------------------------------------------------

Description:
  Detect virulence genes in whole-genome sequences using BLAST search against the FDA VirulenceFinder database (https://virulence.fda.gov/).

Inputs:
  - fasta: Genome assembly in FASTA format
  - species: Bacterial species (e.g. "Nmeningitidis")
  - sample_name: Sample ID or name

Outputs:
  - virulence_hits.csv: CSV file listing matched virulence genes
  - virulence_report.html: HTML table of detected virulence genes
  - logs/: Log folder containing BLAST and preprocessing logs

Command:
  miniwdl run virulence_finder.wdl fasta=assembly.fasta species=Nmeningitidis sample_name=Sample123

Notes:
  - Uses a preformatted JSON file with known virulence profiles per species.
  - Requires BLAST+ and pandas with Jinja2 installed in the environment.

-------------------------------------------------------------------------------
3. Workflow: quality_control.wdl
-------------------------------------------------------------------------------

Description:
  Perform quality control on raw reads using FastQC and generate a summary report.

Inputs:
  - read1: First FASTQ file (paired-end)
  - read2: Second FASTQ file (paired-end)

Outputs:
  - fastqc_read1.html / .zip: FastQC report and data for Read 1
  - fastqc_read2.html / .zip: FastQC report and data for Read 2
  - summary.txt: Aggregated QC summary

Command:
  miniwdl run quality_control.wdl read1=sample_R1.fastq.gz read2=sample_R2.fastq.gz

Notes:
  - FastQC must be accessible in the Docker container or host machine.

-------------------------------------------------------------------------------
4. Workflow: assemble_spades.wdl
-------------------------------------------------------------------------------

Description:
  Assemble paired-end reads using SPAdes assembler.

Inputs:
  - read1: First FASTQ file
  - read2: Second FASTQ file

Outputs:
  - contigs.fasta: Final genome assembly
  - assembly_graph.gfa: Assembly graph
  - spades.log: Log file from SPAdes

Command:
  miniwdl run assemble_spades.wdl read1=sample_R1.fastq.gz read2=sample_R2.fastq.gz

Notes:
  - Default SPAdes options are used; parameters can be tuned by modifying the task.

-------------------------------------------------------------------------------
General Notes:
-------------------------------------------------------------------------------

- All workflows assume Docker is correctly installed and accessible from WDL runtime.
- Input paths must be absolute or relative to the working directory.
- Logs are stored inside the `work/` folder by default and can help in debugging.

Installation Tips:
------------------

To use MiniWDL:
  pip install miniwdl
  miniwdl run your_workflow.wdl input=value ...

To reset config:
  mv ~/.config/miniwdl.cfg ~/.config/miniwdl.cfg.bak
  miniwdl configure  # interactive

To use with Docker:
  Ensure your `miniwdl.cfg` contains:
    docker = true

-------------------------------------------------------------------------------
License:
  Ministry Of Health - Jerusalme

Maintainers:
  David Maimoun

Contact:
  For questions or help, please open an issue or contact david.maimoun@moh.gov.il.


