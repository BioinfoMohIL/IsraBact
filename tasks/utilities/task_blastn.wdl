version 1.0

## task_blastn.wdl
## Tasks pour detecter la presence d'un gene dans des assemblages (blastn)

task blastn_gene {
  meta {
    description: "Cree une base BLAST a partir d'un assemblage, blastn le gene, et dit si le gene est present (seuils identite + couverture)."
  }

  input {
    File assembly                     # assemblage du sample (fasta)
    File gene_fasta                   # sequence du gene (fasta)
    String sample_name
    Float min_identity = 90.0         # % identite minimum
    Float min_coverage = 80.0         # % du gene couvert minimum
    String evalue = "1e-10"            # String pour eviter la notation scientifique (non supportee par le parser WDL 1.0)
    String docker = "ncbi/blast:2.15.0"
    Int cpu = 2
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail

    makeblastdb -in ~{assembly} -dbtype nucl -out db_~{sample_name}

    blastn \
      -query ~{gene_fasta} \
      -db db_~{sample_name} \
      -evalue ~{evalue} \
      -num_threads ~{cpu} \
      -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
      -out ~{sample_name}_blast.tsv

    # Colonnes: 1 qseqid 2 sseqid 3 pident 4 length 5 qlen 6 qstart 7 qend 8 sstart 9 send 10 evalue 11 bitscore
    # On garde le meilleur hit (bitscore) qui passe les seuils identite + couverture
    awk -F'\t' -v s="~{sample_name}" -v mi=~{min_identity} -v mc=~{min_coverage} '
      BEGIN { OFS="\t"; best=-1 }
      {
        cov = ($7 - $6 + 1) / $5 * 100
        if ($3 >= mi && cov >= mc && $11 > best) {
          best = $11
          out = s OFS "PRESENT" OFS $3 OFS sprintf("%.1f", cov) OFS $2 OFS $8 OFS $9 OFS $10 OFS $11
        }
      }
      END {
        if (best >= 0) print out
        else print s, "ABSENT", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
      }' ~{sample_name}_blast.tsv > ~{sample_name}_summary.tsv

    cut -f2 ~{sample_name}_summary.tsv > gene_status.txt
  >>>

  output {
    File blast_tsv = "~{sample_name}_blast.tsv"
    File summary = "~{sample_name}_summary.tsv"
    String gene_status = read_string("gene_status.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
  }
}

task merge_summaries {
  meta {
    description: "Concatene les resumes par sample en un seul tableau avec en-tete."
  }

  input {
    Array[File] summaries
    String output_name = "gene_presence_summary.tsv"
    String docker = "ubuntu:22.04"
  }

  command <<<
    set -euo pipefail
    printf "sample\tgene_status\tpident\tgene_coverage\tcontig\tsstart\tsend\tevalue\tbitscore\n" > ~{output_name}
    cat ~{sep=" " summaries} >> ~{output_name}
  >>>

  output {
    File merged = "~{output_name}"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
  }
}
