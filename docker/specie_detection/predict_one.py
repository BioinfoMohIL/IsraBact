"""Prédit l'espèce d'UN échantillon paired-end (utilisé par la task WDL ml_species_predict).

Écrit dans --out_dir: species.txt, top_taxa.txt, version.txt, report.tsv
"""
import argparse
import hashlib
import pickle
import sys
import warnings
from pathlib import Path

import sklearn

from kmer_utils import sample_features

MIN_READS_WARN = 10_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--sample_id", default="sample")
    ap.add_argument("--threshold", type=float, help="Écrase le seuil enregistré dans le modèle")
    ap.add_argument("--out_dir", default=".")
    args = ap.parse_args()

    model_bytes = Path(args.model).read_bytes()
    b = pickle.loads(model_bytes)
    if b["sklearn_version"] != sklearn.__version__:
        warnings.warn(f"Modèle entraîné avec sklearn {b['sklearn_version']}, version actuelle {sklearn.__version__}")
    model, k, max_reads = b["model"], b["k"], b["max_reads"]
    if len(model.classes_) < 2:
        sys.exit(f"Modèle à une seule classe ({list(model.classes_)}): ré-entraîner avec au moins 2 espèces.")
    threshold = args.threshold if args.threshold is not None else b["threshold"]

    # EXACTEMENT la même featurisation que pendant l'entraînement
    x, n_reads = sample_features(args.r1, args.r2, k=k, max_reads=max_reads)
    if n_reads == 0:
        sys.exit("Aucun read lu dans les FASTQ (fichiers vides ?).")
    proba = model.predict_proba(x.reshape(1, -1))[0]
    order = proba.argsort()[::-1]
    classes = [str(c) for c in model.classes_]

    best, p1 = classes[order[0]], float(proba[order[0]])
    species = best if p1 >= threshold else "UNKNOWN"
    qc = "LOW_READS" if n_reads < MIN_READS_WARN else "ok"
    top_taxa = "; ".join(f"{classes[i]} ({proba[i]:.2f})" for i in order[:3])
    model_id = hashlib.md5(model_bytes).hexdigest()[:8]

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "species.txt").write_text(species)
    (out / "top_taxa.txt").write_text(top_taxa)
    # genre seul, pour task_species_match ("" si UNKNOWN: rien de détecté, comme kraken2 sans classification)
    genus = best.split()[0] if species != "UNKNOWN" else ""
    (out / "genus.txt").write_text(genus + "\n")
    (out / "version.txt").write_text(f"strain_clf RandomForest k={k} model={model_id} sklearn={sklearn.__version__}")
    with open(out / "report.tsv", "w") as f:
        f.write(f"# sample={args.sample_id} n_reads_used={n_reads} qc={qc} threshold={threshold} "
                f"best_guess={best} confidence={p1:.3f}\n")
        f.write("species\tprobability\n")
        for i in order:
            f.write(f"{classes[i]}\t{proba[i]:.4f}\n")
    print(f"{args.sample_id}: {species} (best={best}, p={p1:.2f}) [{qc}, {n_reads} reads]")


if __name__ == "__main__":
    main()
