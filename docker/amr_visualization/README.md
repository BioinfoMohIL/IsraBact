# AMR heatmap — phylogénie + résistome/virulome sur Terra

Recrée la **Fig. 3e** d'AMRViz (heatmap résistance/virulence alignée sur une
phylogénie) à partir de rapports **Abricate** / **ResFinder** déjà produits sur
Terra. **Python pur** (pandas, matplotlib, biopython) + MAFFT/IQ-TREE pour
l'arbre. Aucune web app, aucun port ouvert → s'exécute proprement en WDL.

---

## Pourquoi pas AMRViz directement sur Terra ? (mon avis)

AMRViz = **AMRomics** (pipeline Python) **+ web-app Vue/NodeJS** qui sert un
dashboard sur le **port 3000** (`amrviz.py start`). Les deux questions sont à
séparer :

**1. Le pipeline d'analyse** → oui, faisable sur Terra. C'est de l'assemblage,
annotation (Prokka), BLAST contre les bases AMR/VF, pangénome (panta), phylogénie
(IQ-TREE). Tout ça est du compute batch, parfaitement encapsulable en WDL.

**2. Le dashboard web (port 3000)** → là est le vrai blocage, et il est plus
nuancé que « impossible ». Terra (via Leonardo) **sait** proxy-fier une app web
tournant dans un environnement interactif Jupyter/RStudio — donc techniquement on
*peut* exposer un port. **Mais** pour AMRViz précisément c'est fragile et je ne
le recommande pas :

- l'app Vue/Node a des chemins/base-URL en dur ; passer derrière le proxy
  Leonardo (réécriture d'URL, websockets) casse en général le routage ;
- WDL/Cromwell (le cœur batch de Terra) **n'expose aucun port** — le proxy
  n'existe que dans les VM interactives, pas dans les workflows ;
- ce serait un montage non supporté, à reconstruire à chaque MAJ.

**Verdict** : faire tourner *tout* AMRViz sur Terra n'a pas de bon ROI. L'archi
saine est celle qu'AMRViz permet déjà par construction : **compute découplé de la
visualisation**. Et comme tu as **déjà** les sorties en amont (assemblies,
Abricate, ResFinder), inutile de rejouer le pipeline : on régénère juste la
figure. C'est ce que fait ce package.

> Si un jour tu veux quand même le dashboard interactif : exporte les
> `web-app/static/data/*` vers GCS et sers l'app **ailleurs** (un VM perso, un
> Cloud Run, ou en local si l'IT le permettait). Pas dans un WDL.

---

## Image Docker (Docker Hub)

Image **sans conda/mamba** : `python:3.11-slim` + MAFFT (apt) + binaire statique
IQ-TREE2. Build sur ta machine Linux, push sur Docker Hub.

> **Important — lance `docker build` depuis la racine du repo** (le dossier qui
> contient à la fois `Dockerfile` et `scripts/`). L'erreur `"/scripts": not
> found` vient de là : si le contexte de build ne contient pas `scripts/`, le
> `COPY scripts/ ...` échoue.

```bash
cd amr-heatmap                         # racine : contient Dockerfile + scripts/
docker build -t bioinfomoh/amr-heatmap:latest .
# sur Mac/ARM uniquement, force l'archi amd64 attendue par Terra :
# docker build --platform=linux/amd64 -t bioinfomoh/amr-heatmap:latest .

docker login                           # compte bioinfomoh
docker push bioinfomoh/amr-heatmap:latest
```

> Le repo Docker Hub doit être **public** pour que Terra le tire sans
> authentification. Sinon, il faut configurer les credentials du registre côté
> Terra (plus lourd).

Transfert des fichiers : Terra lit/écrit directement sur GCS (`gs://`), donc pas
de `gsutil` en local.

---

## Architecture

```
            tasks/task_amr_visualization.wdl        workflows/
  ┌─────────────────────────────────────┐   ┌──────────────────────────┐
  │  task ParseAmr                       │   │ amr_heatmap.wdl          │
  │  task BuildTree   (optionnel)        │←──│  import ../tasks/        │
  │  task PlotHeatmap                    │   │   task_amr_visualization │
  └─────────────────────────────────────┘   │  ParseAmr → BuildTree?   │
                                             │          → PlotHeatmap   │
  scripts/ : parse_amr.py · plot_heatmap.py  └──────────────────────────┘
             (logique métier, testable seule)
```

Les trois **tasks** sont regroupées dans `tasks/task_amr_visualization.wdl` ; le
**workflow** ne fait qu'importer ce fichier (`as tasks`) et câbler les appels
(`tasks.ParseAmr`, `tasks.BuildTree`, `tasks.PlotHeatmap`).

---

## Flux WDL

1. **ParseAmr** — Abricate + ResFinder → `matrix.tsv` (échantillons × gènes) +
   `genes.tsv` (catégorie résistance/virulence par gène).
2. **BuildTree** *(conditionnelle)* — si pas d'arbre fourni mais un `core_fasta`
   donné : MAFFT puis IQ-TREE2 → `.nwk`. Si tu as déjà l'arbre, cette task est
   sautée.
3. **PlotHeatmap** — figure phylogénie + heatmap (PDF + PNG), bandeau de
   catégorie de gènes, annotations latérales (ST, espèce…).

---

## Utilisation sur Terra

1. Build + push l'image sur Docker Hub (voir section ci-dessus).
2. **Dockstore** : connecte le repo, publie le workflow `workflows/amr_heatmap.wdl`
   (Dockstore conserve l'arborescence, donc l'import relatif
   `../tasks/task_amr_visualization.wdl` fonctionne), puis « Export to Terra ».
3. Dans Terra, renseigne les inputs (voir `inputs.example.json`) : `docker` =
   `bioinfomoh/amr-heatmap:latest`, pointe tes `gs://…` Abricate/ResFinder, et
   soit `tree` (arbre prêt) soit `core_fasta` (pour lancer IQ-TREE).
4. Lance → Terra écrit `*_heatmap.pdf` / `.png` dans le bucket du workspace.

---

## Lancer en local (test / debug)

```bash
# 1) matrice
python scripts/parse_amr.py \
  --abricate test_data/abricate_resfinder.tsv test_data/vfdb/*.tsv \
  --resfinder test_data/resfinder/*.tsv \
  --min-identity 90 --min-coverage 80 --mode presence \
  --out test_data/matrix.tsv

# 2) figure
python scripts/plot_heatmap.py \
  --matrix test_data/matrix.tsv --genes test_data/matrix.tsv.genes.tsv \
  --tree test_data/tree.nwk \
  --metadata test_data/metadata.tsv --annotate ST species \
  --out amr_heatmap.pdf --png amr_heatmap.png
```

---

## Fichiers

| Fichier | Rôle |
|---|---|
| `scripts/parse_amr.py` | Abricate + ResFinder → matrice + annotation gènes |
| `scripts/plot_heatmap.py` | Phylogénie + heatmap (style Fig 3e) |
| `tasks/task_amr_visualization.wdl` | les 3 tasks : ParseAmr, BuildTree, PlotHeatmap |
| `workflows/amr_heatmap.wdl` | workflow (import + câblage) |
| `Dockerfile` | image (python:slim + mafft via apt + binaire iqtree2) |
| `inputs.example.json` | gabarit d'inputs Terra |

---

## Notes / limites

- **Appariement des IDs** : l'ID d'échantillon vient de la colonne `#FILE`
  d'Abricate et du nom de fichier ResFinder, normalisés (suffixes `.fasta`,
  `_contigs`, `_resfinder`… retirés). Les feuilles de l'arbre doivent matcher ces
  IDs ; le script signale les non-appariés.
- **`mode identity`** colore la heatmap par % d'identité (viridis) au lieu de
  présence/absence binaire.
- **Core genome** : `BuildTree` aligne le FASTA fourni tel quel. Pour un vrai
  core-genome il faut un pangénome (panta/roary) en amont — hors périmètre ici.
- **VFDB** géré : les gènes de virulence (db `vfdb`) sont taggés `virulence` et
  apparaissent dans un bloc séparé du bandeau de catégorie.
