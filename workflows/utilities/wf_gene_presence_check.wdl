version 1.0

## wf_gene_presence_check.wdl
## Verifie si un gene (sequence nucleotidique) est present dans UN assemblage bacterien.
## Generique : fonctionne pour n'importe quel gene / organisme.
## Concu pour Terra : un run par sample (data table), en parallele.

import "../../tasks/utilities/task_blastn.wdl" as blast

workflow wf_gene_presence_check {
  meta {
    description: "BLASTn alignment of a gene against a sample assembly, returning a PRESENT/ABSENT result based on identity and coverage thresholds."
  }

  input {
    File assembly                     # assemblage du sample (fasta, non compresse)
    File gene_fasta                   # sequence nucleotidique du gene a chercher
    String? sample_name               # optionnel, sinon deduit du nom du fichier
    Float min_identity = 90.0         # % identite minimum (a ajuster selon le gene/la question)
    Float min_coverage = 80.0         # % du gene couvert minimum
  }

  String name = select_first([sample_name, sub(basename(assembly), "\\.(fasta|fa|fna|fas)$", "")])

  call blast.blastn_gene {
    input:
      assembly = assembly,
      gene_fasta = gene_fasta,
      sample_name = name,
      min_identity = min_identity,
      min_coverage = min_coverage
  }

  output {
    String  gene_pres_gene_status = blastn_gene.gene_status
    File    gene_pres_summary_tsv = blastn_gene.summary
    File    gene_pres_blast_tsv = blastn_gene.blast_tsv
  }
}
