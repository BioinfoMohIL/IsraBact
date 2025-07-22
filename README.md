IsraBact WDL Workflows (v1.0)
=======================================

De-novo genome assembly, taxonomic ID, genes and species typing, alleles comparison, and virulence prediction from paired-end bacterial NGS data.

This project contains WDL (Workflow Description Language) workflows based on Theiagen - Public Health Bioinformatics pipelines and tools such as Kraken2, BLAST, and ListPred.

These WDLs are designed to be executed via MiniWDL or Cromwell with Docker backend support.

-------------------------------------------------------------------------------
1. Workflow: metagenomics/wf_species_detection_bs_reads.wdl
-------------------------------------------------------------------------------

Description:
  Check if the sequenced samples on Basespace Illumina correspond to the expected species using Kraken2 classification.

Inputs:
  - String: `api_server`  
      → Illumina Basespace API endpoint  
  - String: `access_token`  
      → Your personal token to access Basespace data  
  - String: `basespace_collection_id`  
      → Name or ID of the sample (e.g., EC001, NM005)  
  - String: `sample_prefix`  
      → Optional prefix to filter samples (e.g., 'EC' for E. coli)

Outputs:
  - File: `reads_list.txt`  
      → List of input reads fetched from Basespace  
  - Array[String]: `samples`  
      → Sample names included in the analysis  
  - File: `species_detected.csv`  
      → Table with sample name, Kraken2-detected species, and match status ("+" if detected species matches prefix-based expected species)


-------------------------------------------------------------------------------
2. Workflow: virulence/wf_virulence_finder.wdl
-------------------------------------------------------------------------------

Description:
  Detects virulence genes using BLAST against the FDA VirulenceFinder database (https://virulence.fda.gov/).  
  Generates CSV and HTML reports showing presence/absence of genes.

Inputs:
  - File: `fasta`  
      → Genome assembly in FASTA format  
  - String: `species`  
      → 'salmonella' or 'ecoli'  
  - String: `samplename`  
      → Sample ID (used in output naming)  
  - Boolean: `plasmid_check` (default: false)  
      → If true, will also search for plasmid virulence genes

Outputs:
  - File: `virfind_hits_csv`  
      → Table with detailed BLAST hits  
  - File: `virfind_summary_csv`  
      → Matrix of gene presence/absence (0/1)  
  - File: `virfind_summary_html`  
      → Visual HTML summary table  
  - Optional Files (if `plasmid_check=true`):  
      → `virfind_plasmid_hits_csv`  
      → `virfind_plasmid_summary_csv`  
      → `virfind_plasmid_summary_html`

Notes:
  - All searches use local BLAST with JSON virulence references.


-------------------------------------------------------------------------------
3. Workflow: resistance/wf_listeria_pred_fasta.wdl
-------------------------------------------------------------------------------

Description:
  Predicts *Listeria monocytogenes* serotypes and phenotypes using ListPred from assembled genomes.

Inputs:
  - File: `assembly_fasta`  
      → Assembled genome in FASTA format  
  - Int: `cpu` (default = 12)  
      → Number of CPU threads for Snakemake

Outputs:
  - File: `virulence_prediction_out.csv`  
  - File: `combined_predictions_out_categorical.csv`  
  - File: `combined_predictions_out_numerical.csv`  
  - File: `disinftolerance_prediction_out.csv`  
  - File: `vir_align_out.tar.gz` (optional)  
  - File: `disinf_align_out.tar.gz` (optional)  
  - String: `virulence_class`  
  - String: `disinfectant_phenotype`

Notes:
  - Outputs are extracted from prediction and alignment results


-------------------------------------------------------------------------------
4. Workflow: resistance/wf_listeria_pred_reads.wdl
-------------------------------------------------------------------------------

Description:
  Predicts *Listeria monocytogenes* serotypes and phenotypes directly from paired-end raw FASTQ reads using ListPred.

Inputs:
  - File: `read1`  
  - File: `read2`  
  - Int: `cpu` (default = 20)

Outputs:
  - File: `virulence_prediction_out.csv`  
  - File: `combined_predictions_out_categorical.csv`  
  - File: `combined_predictions_out_numerical.csv`  
  - File: `disinftolerance_prediction_out.csv`  
  - File: `vir_align_out.tar.gz` (optional)  
  - File: `disinf_align_out.tar.gz` (optional)  
  - String: `virulence_class`  
  - String: `disinfectant_phenotype`

Notes:
  - Input FASTQ files are renamed and reformatted automatically
  - Snakemake handles serotype/virulence/disinfectant classification


-------------------------------------------------------------------------------
5. Workflow: utilities/wf_fq2dna
-------------------------------------------------------------------------------

Description:
  De-novo genome assembly and quality control of paired-end NGS reads using fq2DNA  
  (Developed by Alexis Criscuolo, Institut Pasteur). This workflow processes reads to produce assembly fasta files, metrics, and scaffold/contig selections.

Inputs:
  - File: `read1`  
      → Forward FASTQ reads (paired-end)  
  - File: `read2`  
      → Reverse FASTQ reads (paired-end)  
  - String: `species`  
      → Species identifier used for annotation or filtering  
  - String: `organism` (default: "B")  
      → Organism type code:  
        - "B" = Bacteria  
        - "P" = Prokaryote  
        - "E" = Eukaryote  
        - "S" = Standard  
        - "V" = Virus  
  - String: `alien_tag` (default: "AUTO")  
      → Tag for contamination or foreign DNA detection

Outputs:
  - String: `fq2dna_version`  
      → fq2DNA software version used  
  - File: `fq2dna_assembly_fasta`  
      → Genome assembly FASTA file produced by fq2DNA  
  - File: `fq2dna_metrics_zip`  
      → ZIP archive containing quality metrics and reports  
  - File: `fq2dna_selected_scaffolds`  
      → Selected scaffold sequences in FASTA format  
  - File: `fq2dna_selected_contigs`  
      → Selected contigs in FASTA format  
  - File: `fq2dna_scaffolding_info`  
      → Text file summarizing scaffolding information and statistics

Notes:
  - This workflow relies on the `fq2dna_run` task which executes the fq2DNA pipeline.
  - Organism and alien_tag parameters help tune the assembly and contamination filtering.
  - Outputs include both assembled sequences and detailed QC metrics for downstream analysis.


-------------------------------------------------------------------------------
6. Workflow: utilities/assembly_long_reads
-------------------------------------------------------------------------------

Description:
  De-novo genome assembly for long reads primarily using Unicycler with optional Pilon polishing,
  filtering of contigs, and quality control metrics generation (QUAST, BUSCO).
  Includes trimming and QC of paired-end reads before assembly.

Inputs:
  - File: `read1`  
      → Forward reads FASTQ (paired-end)  
  - File: `read2`  
      → Reverse reads FASTQ (paired-end)  
  - String: `samplename`  
      → Sample identifier  
  - String: `assembler` (default: "unicycler")  
      → Assembly tool to use (currently supports Unicycler)  
  - String: `seq_method` (default: "Illumina")  
      → Sequencing method description  
  - Int: `trim_min_length` (default: 75)  
      → Minimum read length after trimming  
  - Int: `trim_quality_min_score` (default: 20)  
      → Minimum base quality score for trimming  
  - Int: `trim_window_size` (default: 10)  
      → Sliding window size for trimming  
  - File?: `illumina_unpaired_fq`  
      → Optional unpaired Illumina reads  
  - File?: `long_reads`  
      → Optional long reads (e.g., Nanopore, PacBio) for hybrid assembly  
  - Boolean: `call_pilon` (default: false)  
      → Whether to run Pilon polishing post-assembly  
  - String?: `assembler_options`  
      → Extra assembler options as a string  
  - Boolean: `run_filter_contigs` (default: false)  
      → Whether to filter contigs after assembly  
  - Additional optional parameters for bwa, pilon, and filtering tasks (CPU, memory, docker images, etc.)

Outputs:
  - Multiple outputs related to:
    - Read QC (FastQC, FastQ screen, trimming reports, cleaned reads)
    - Assembly results (FASTA, GFA, logs)
    - Pilon polishing outputs (VCF, changes)
    - Filtered contigs metrics
    - Assembly QC reports (QUAST, BUSCO), including assembly length, N50, GC content, completeness scores

Notes:
  - Performs paired-end read trimming and quality control before assembly.
  - Supports hybrid assembly with optional long reads input.
  - Pilon polishing improves assembly accuracy if enabled.
  - Filtering step removes low-quality or short contigs.
  - Generates detailed QC and assembly reports for downstream analysis.


-------------------------------------------------------------------------------

Environment:
  - Requires Docker
  - MiniWDL or Cromwell recommended
  - Internet access required for downloading containers or accessing Basespace

License:
  Ministry of Health – Jerusalem, Israel

Version:
  1.0

Maintainer:
  David Maimoun  
  Email: david.maimoun@moh.gov.il

Support:
  Please open an issue on the repository or email the maintainer directly.
