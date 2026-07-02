# Données d'exemple — 8 isolats *K. pneumoniae* MDR (synthétiques)

Jeu de données factice pour tester le pipeline de bout en bout.

```
bash examples/run_example.sh
# -> examples/outputs/KpMDR_heatmap.pdf + .png
```

## Format des fichiers d'entrée

### `inputs/abricate_resfinder.tsv` — résistome (Abricate combiné)
Sortie Abricate standard, **plusieurs échantillons** dans un seul fichier,
distingués par la colonne `#FILE`. L'ID échantillon est dérivé de `#FILE`
(`Kp01.fasta` → `Kp01`).

| #FILE | SEQUENCE | START | END | STRAND | GENE | %COVERAGE | %IDENTITY | DATABASE | ACCESSION | PRODUCT |
|---|---|---|---|---|---|---|---|---|---|---|

### `inputs/vfdb/*.tsv` — virulome (Abricate, base VFDB)
Même format Abricate, **un fichier par échantillon**, `DATABASE = vfdb`. Ces
gènes sont automatiquement taggés `virulence` et placés dans un bloc séparé.

### `inputs/resfinder/*.tsv` — ResFinder (un fichier par échantillon)
Sortie `ResFinder_results_tab.txt`. L'ID vient du nom de fichier
(`Kp01_resfinder.tsv` → `Kp01`).

| Resistance gene | Identity | Alignment Length/Gene Length | Coverage | ... | Phenotype | Accession no. |
|---|---|---|---|---|---|---|

### `inputs/tree.nwk` — phylogénie (Newick)
Arbre pré-calculé. **Les noms de feuilles doivent matcher les IDs** (`Kp01`…).
Si tu n'as pas d'arbre, fournis plutôt un `core_fasta` au workflow → IQ-TREE2.

### `inputs/metadata.tsv` — annotations (optionnel)
Colonne `sample` + colonnes à afficher à côté de l'arbre (ici `ST`, `species`).

| sample | ST | species |
|---|---|---|
| Kp01 | ST258 | K.pneumoniae |

## Sorties (`outputs/`)
- `KpMDR_matrix.tsv` — matrice échantillons × gènes (présence/absence)
- `KpMDR_matrix.tsv.genes.tsv` — catégorie par gène (resistance/virulence)
- `KpMDR_heatmap.pdf` / `.png` — figure phylogénie + heatmap

> Données **synthétiques** (générées aléatoirement) : la topologie de l'arbre et
> les profils de gènes n'ont aucune signification biologique, ils servent
> seulement à valider le format et le rendu.
