#!/usr/bin/env python3
"""
parse_amr.py - Build a gene x sample matrix from resistance and (optional)
virulence reports.

Resistance input (mandatory): one file per sample, in EITHER format - the
script auto-detects each file:
  - AMRFinderPlus  (columns 'Name', 'Element symbol', 'Type', ...)
  - ResFinder      ('ResFinder phenotype results.' header, OR a results_tab)

Virulence input (optional): Abricate reports run against a virulence DB (vfdb).

Outputs:
  - <out>            : TSV matrix (rows = samples, cols = genes),
                       values = 1/0 (presence) or %identity (--mode identity)
  - <out>.genes.tsv  : column annotation (gene, category, drug_class)

The matrix and annotation feed plot_heatmap.py directly, which then has no
parsing logic of its own.
"""
import argparse
import os
import sys
import pandas as pd

# Common assembly/file suffixes stripped to derive a clean sample id.
STRIP_SUFFIXES = [
    "_resfinder", ".resfinder", "_abricate", ".abricate", "_amrfinder",
    "_amr", "_all", "_hits", "_results_tab", "_pheno_table", "_phenotype",
    ".tsv", ".txt", ".tab", ".tab.txt",
    "_filtered_contigs", "_filtered", "_contigs", "_assembly", "_spades",
    "_skesa", "_shovill", "_polished",
    ".fasta", ".fa", ".fna", ".contigs", ".assembly",
]


def normalize_id(raw: str) -> str:
    """Derive a clean sample id from a filename / path / FILE field."""
    name = os.path.basename(str(raw)).strip()
    changed = True
    while changed:
        changed = False
        low = name.lower()
        for suf in STRIP_SUFFIXES:
            if low.endswith(suf):
                name = name[: len(name) - len(suf)]
                changed = True
                low = name.lower()
    return name or os.path.basename(str(raw))


# --- format detection --------------------------------------------------------
def detect_format(path):
    """Return one of: amrfinder, resfinder_pheno, resfinder_tab, abricate, unknown."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [fh.readline() for _ in range(40)]
    nonempty = [l for l in lines if l.strip()]
    if not nonempty:
        return "unknown"
    if nonempty[0].lower().startswith("# resfinder phenotype"):
        return "resfinder_pheno"
    for l in nonempty:
        low = l.lower()
        if low.startswith("#file\t") or low.startswith("#file "):
            return "abricate"
        if not l.startswith("#"):
            cols = [c.strip().lower() for c in l.rstrip("\n").split("\t")]
            if "element symbol" in cols:
                return "amrfinder"
            if "resistance gene" in cols:
                return "resfinder_tab"
            if "gene" in cols and "database" in cols:
                return "abricate"
            break
    return "unknown"


# --- individual parsers ------------------------------------------------------
# Each returns (records, samples_seen).
# record = (sample, gene, identity, category, drug_class)
AMRFINDER_CATEGORY = {"AMR": "resistance", "VIRULENCE": "virulence", "STRESS": "stress"}


def parse_amrfinder(path, min_id, min_cov):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cols = {c.lower().strip(): c for c in df.columns}
    name_c = cols.get("name")
    gene_c = cols.get("element symbol")
    id_c = cols.get("% identity to reference")
    cov_c = cols.get("% coverage of reference")
    type_c = cols.get("type")
    class_c = cols.get("class")

    samples, records = set(), []
    for _, row in df.iterrows():
        sample = normalize_id(row[name_c]) if name_c else normalize_id(path)
        samples.add(sample)
        ident = float(row[id_c]) if id_c and row[id_c] not in ("", "NA") else 100.0
        cov = float(row[cov_c]) if cov_c and row[cov_c] not in ("", "NA") else 100.0
        if ident < min_id or cov < min_cov:
            continue
        gene = str(row[gene_c]).strip()
        if not gene:
            continue
        typ = str(row[type_c]).strip().upper() if type_c else "AMR"
        category = AMRFINDER_CATEGORY.get(typ, "resistance")
        drug_class = str(row[class_c]).strip() if class_c else ""
        records.append((sample, gene, ident, category, drug_class))
    return records, samples


def parse_resfinder_pheno(path, min_id, min_cov):
    """ResFinder phenotype table: genes live in the 'Genetic background' column."""
    sample = normalize_id(path)
    records = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("# Sample:"):
                sample = normalize_id(line.split(":", 1)[1])
            if line.startswith("#") or not line.strip():
                continue
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            if len(parts) < 5:
                continue
            drug_class, genetic_bg = parts[1], parts[4]
            if not genetic_bg:
                continue
            # tokens like "blaTEM-1B (blaTEM-1B_AY458016)", comma-separated
            for tok in genetic_bg.split(","):
                gene = tok.strip().split(" (")[0].strip()
                if gene:
                    records.append((sample, gene, 100.0, "resistance", drug_class))
    return records, {sample}


def parse_resfinder_tab(path, min_id, min_cov):
    """Older ResFinder 'ResFinder_results_tab.txt' (one gene per row)."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cols = {c.lower().strip(): c for c in df.columns}
    gene_c = cols.get("resistance gene") or cols.get("gene")
    id_c = cols.get("identity")
    cov_c = cols.get("coverage")
    class_c = cols.get("phenotype")
    sample = normalize_id(path)
    records = []
    for _, row in df.iterrows():
        ident = float(row[id_c]) if id_c and row[id_c] not in ("", "NA") else 100.0
        cov = float(row[cov_c]) if cov_c and row[cov_c] not in ("", "NA") else 100.0
        if ident < min_id or cov < min_cov:
            continue
        gene = str(row[gene_c]).strip()
        if gene:
            records.append((sample, gene, ident, "resistance",
                            str(row[class_c]).strip() if class_c else ""))
    return records, {sample}


def parse_abricate(path, min_id, min_cov, default_category="virulence"):
    """Abricate report. The #FILE column is the authoritative sample id."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df.columns = [c.strip().lstrip("#").upper() for c in df.columns]
    if "GENE" not in df.columns:
        sys.stderr.write(f"[abricate] no GENE column in {path}, skipped\n")
        return [], set()
    if "FILE" in df.columns:
        df["__sample"] = df["FILE"].map(normalize_id)
    else:
        df["__sample"] = normalize_id(path)

    pid = pd.to_numeric(df.get("%IDENTITY", 100), errors="coerce").fillna(100.0)
    pcov = pd.to_numeric(df.get("%COVERAGE", 100), errors="coerce").fillna(100.0)
    db_cat = {"vfdb": "virulence", "ecoli_vf": "virulence",
              "resfinder": "resistance", "card": "resistance",
              "ncbi": "resistance", "argannot": "resistance",
              "plasmidfinder": "plasmid"}

    samples, records = set(), []
    for sample, gene, ident, cov, db in zip(
            df["__sample"], df["GENE"], pid, pcov,
            df.get("DATABASE", pd.Series([""] * len(df)))):
        samples.add(sample)
        if ident < min_id or cov < min_cov:
            continue
        category = db_cat.get(str(db).lower().strip(), default_category)
        records.append((sample, str(gene).strip(), float(ident), category, ""))
    return records, samples


# --- dispatch ----------------------------------------------------------------
def parse_resistance_file(path, min_id, min_cov):
    fmt = detect_format(path)
    if fmt == "amrfinder":
        return parse_amrfinder(path, min_id, min_cov)
    if fmt == "resfinder_pheno":
        return parse_resfinder_pheno(path, min_id, min_cov)
    if fmt == "resfinder_tab":
        return parse_resfinder_tab(path, min_id, min_cov)
    if fmt == "abricate":  # abricate run against a resistance DB
        return parse_abricate(path, min_id, min_cov, default_category="resistance")
    sys.stderr.write(f"[warn] unrecognized resistance format: {path} (skipped)\n")
    return [], set()


def parse_virulence_file(path, min_id, min_cov):
    fmt = detect_format(path)
    if fmt == "abricate":
        return parse_abricate(path, min_id, min_cov, default_category="virulence")
    if fmt == "amrfinder":  # AMRFinder '--plus' may carry VIRULENCE rows
        recs, smp = parse_amrfinder(path, min_id, min_cov)
        return [r for r in recs if r[3] == "virulence"], smp
    sys.stderr.write(f"[warn] unrecognized virulence format: {path} (skipped)\n")
    return [], set()


# --- matrix ------------------------------------------------------------------
def build_matrix(records, all_samples, mode="presence"):
    df = pd.DataFrame(records, columns=["sample", "gene", "identity", "category", "drug_class"])
    if df.empty:
        matrix = pd.DataFrame(index=sorted(all_samples))
        ann = pd.DataFrame(columns=["gene", "category", "drug_class"])
        return matrix, ann

    df = df.sort_values("identity").drop_duplicates(["sample", "gene"], keep="last")
    df["val"] = 1.0 if mode == "presence" else df["identity"]
    matrix = df.pivot_table(index="sample", columns="gene", values="val", fill_value=0.0)

    # include samples that had zero genes (e.g. a 'No resistance' ResFinder table)
    matrix = matrix.reindex(sorted(set(matrix.index) | set(all_samples)), fill_value=0.0)

    ann = (df.groupby("gene")
             .agg(category=("category", lambda s: s.value_counts().index[0]),
                  drug_class=("drug_class", lambda s: next((x for x in s if x), "")))
             .reindex(matrix.columns).reset_index())
    return matrix, ann


def main():
    ap = argparse.ArgumentParser(description="Resistance/virulence reports -> gene x sample matrix")
    ap.add_argument("--resistance", nargs="+", required=True,
                    help="resistance reports (AMRFinder or ResFinder; one per sample)")
    ap.add_argument("--virulence", nargs="*", default=[],
                    help="virulence reports (Abricate vfdb); optional")
    ap.add_argument("--min-identity", type=float, default=90.0)
    ap.add_argument("--min-coverage", type=float, default=80.0)
    ap.add_argument("--mode", choices=["presence", "identity"], default="presence")
    ap.add_argument("--out", required=True, help="output matrix path (TSV)")
    args = ap.parse_args()

    records, all_samples = [], set()
    for p in args.resistance:
        recs, smp = parse_resistance_file(p, args.min_identity, args.min_coverage)
        records += recs
        all_samples |= smp
    for p in args.virulence:
        recs, smp = parse_virulence_file(p, args.min_identity, args.min_coverage)
        records += recs
        all_samples |= smp

    if not all_samples:
        sys.exit("No samples parsed: check the resistance inputs and thresholds.")

    matrix, genes_ann = build_matrix(records, all_samples, args.mode)
    matrix.index.name = "sample"
    matrix.to_csv(args.out, sep="\t")
    genes_ann.to_csv(args.out + ".genes.tsv", sep="\t", index=False)

    n_res = sum(1 for r in records if r[3] == "resistance")
    n_vir = sum(1 for r in records if r[3] == "virulence")
    sys.stderr.write(
        f"[ok] {matrix.shape[0]} samples x {matrix.shape[1]} genes -> {args.out}\n"
        f"     resistance hits: {n_res} | virulence hits: {n_vir} | "
        f"other: {len(records) - n_res - n_vir}\n"
        f"     annotation -> {args.out}.genes.tsv\n")


if __name__ == "__main__":
    main()
