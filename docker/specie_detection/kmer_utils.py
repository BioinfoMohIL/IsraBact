"""Fonctions partagées entre l'entraînement et la prédiction.

Une seule source de vérité pour la featurisation = pas de "train/predict skew"
(le piège classique où le modèle est entraîné sur des features calculées d'une
façon et testé sur des features calculées d'une autre).
"""
import gzip
import re
from functools import lru_cache
from itertools import islice
from pathlib import Path

import numpy as np

FASTQ_EXT = (".fastq", ".fastq.gz", ".fq", ".fq.gz")

# Table de correspondance base -> entier. 4 = base inconnue (N, etc.)
_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i


def read_fastq(path, max_reads):
    """Retourne les `max_reads` premières séquences (1 ligne sur 4)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line.strip() for line in islice(f, 1, 4 * max_reads, 4)]


def find_pairs(folder):
    """Trouve les paires R1/R2 dans un dossier.

    Gère: sample_R1_001.fastq.gz / sample_R2_001.fastq.gz  et  sample_1.fastq.gz / sample_2.fastq.gz
    Retourne une liste de (sample_id, chemin_R1, chemin_R2).
    """
    pairs = []
    for p in sorted(Path(folder).iterdir()):
        if not p.name.endswith(FASTQ_EXT):
            continue
        m = re.search(r"_R1(?=[_.])", p.name) or re.search(r"_1(?=\.)", p.name)
        if m is None:
            continue  # ce n'est pas un R1
        r2_name = p.name[:m.start()] + m.group().replace("1", "2") + p.name[m.end():]
        r2 = p.with_name(r2_name)
        if r2.exists():
            pairs.append((p.name[:m.start()], p, r2))
        else:
            print(f"[skip] {p.name}: pas de R2 correspondant")
    return pairs


@lru_cache(maxsize=None)
def canonical_indices(k):
    """Indices (parmi les 4**k k-mers possibles) qui sont 'canoniques'.

    Un k-mer et son reverse-complement représentent le même brin d'ADN vu de
    l'autre côté. On ne garde que le plus petit des deux (en code entier).
    """
    idx = np.arange(4 ** k, dtype=np.int64)
    rc = np.zeros_like(idx)
    for i in range(k):
        digit = (idx // 4 ** (k - 1 - i)) % 4
        rc += (3 - digit) * 4 ** i
    return idx[idx <= rc]


def sample_features(r1, r2, k=8, max_reads=100_000):
    """Vecteur de fréquences de k-mers canoniques pour un échantillon paired-end.

    Retourne (vecteur float32 normalisé qui somme à 1, nombre de reads réellement lus).
    """
    n_feat = len(canonical_indices(k))
    reads = read_fastq(r1, max_reads // 2) + read_fastq(r2, max_reads // 2)
    n_reads = len(reads)
    if n_reads == 0:
        return np.zeros(n_feat, dtype=np.float32), 0

    # Tous les reads collés avec des "N" comme séparateurs: un seul passage numpy
    # (beaucoup plus rapide que Counter sur des strings), et les k-mers qui
    # chevauchent deux reads contiennent un N donc sont ignorés.
    seq = "N".join(reads).upper()
    arr = _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)].astype(np.int64)
    n = len(arr) - k + 1
    if n <= 0:
        return np.zeros(n_feat, dtype=np.float32), n_reads

    fwd = np.zeros(n, dtype=np.int64)    # code du k-mer lu dans le sens normal
    rc = np.zeros(n, dtype=np.int64)     # code de son reverse-complement
    bad = np.zeros(n, dtype=bool)        # k-mer contenant un N
    for i in range(k):
        seg = arr[i:i + n]
        bad |= seg > 3
        fwd += seg * 4 ** (k - 1 - i)
        rc += (3 - seg) * 4 ** i
    canon = np.minimum(fwd, rc)[~bad]

    counts = np.bincount(canon, minlength=4 ** k)[canonical_indices(k)].astype(np.float32)
    total = counts.sum()
    return (counts / total if total > 0 else counts), n_reads
