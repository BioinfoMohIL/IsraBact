# Lancer le workflow en local via `inputs.json`

Tu testes le **workflow WDL complet** sur ta machine Linux (Docker dispo), sans
Terra. Deux runners possibles : `miniwdl` (le plus simple) ou Cromwell.

## 0. Pré-requis (une fois)

```bash
# image construite EN LOCAL -> le runner la trouve sans push Docker Hub
cd amr-heatmap
docker build -t bioinfomoh/amr-heatmap:latest .

pip install miniwdl        # runner leger
```

## Quel inputs choisir ?

Le seul input **obligatoire** est `docker`. Pour la phylogénie, tu fournis
**soit** un arbre **soit** un FASTA — d'où les deux fichiers :

| Fichier | Phylogénie | Ce qui tourne | Quand l'utiliser |
|---|---|---|---|
| `examples/inputs.local.tree.json` | `tree` fourni | ParseAmr → PlotHeatmap (IQ-TREE **sauté**) | **1er test** : rapide, peu de pièces mobiles |
| `examples/inputs.local.fasta.json` | `core_fasta` | ParseAmr → **BuildTree (MAFFT+IQ-TREE)** → PlotHeatmap | tester la branche phylogénie |

Règles de sélection :
- **`tree` OU `core_fasta`**, pas les deux. Si les deux sont donnés, `tree`
  gagne et `BuildTree` ne tourne pas.
- Si **aucun** des deux : le workflow tourne quand même, mais la heatmap est
  triée alphabétiquement (pas de panneau arbre).
- Il faut **au moins** `abricate_reports` ou `resfinder_reports` non vide,
  sinon `ParseAmr` s'arrête (aucun gène).

## 1. Premier test — arbre fourni (recommandé)

Depuis la **racine du repo** (les chemins du JSON sont relatifs au CWD) :

```bash
miniwdl run workflows/amr_heatmap.wdl \
  -i examples/inputs.local.tree.json \
  -d runs/
```

À la fin, miniwdl affiche un JSON des outputs avec les chemins ; la figure est
sous `runs/<...>/out/heatmap_pdf/KpMDR_heatmap.pdf` (et `_heatmap.png`).

## 2. Test branche IQ-TREE — core_fasta

```bash
miniwdl run workflows/amr_heatmap.wdl \
  -i examples/inputs.local.fasta.json \
  -d runs/
```

Ici `BuildTree` aligne `examples/inputs/core.fasta` (MAFFT) puis infère l'arbre
(IQ-TREE2) avant la heatmap. Plus long (quelques minutes).

## Variante Cromwell

```bash
# telecharge cromwell.jar depuis github.com/broadinstitute/cromwell/releases
java -jar cromwell.jar run workflows/amr_heatmap.wdl \
  -i examples/inputs.local.tree.json
```

Les imports relatifs (`../tasks/task_amr_visualization.wdl`) sont résolus par
rapport à l'emplacement du `.wdl` — pas besoin de config supplémentaire.

## Surcharger un paramètre

Édite le JSON, ou passe `-i` avec un second fichier qui surcharge. Exemples de
clés utiles : `AmrHeatmap.mode` = `"identity"` (heatmap colorée par %identité),
`AmrHeatmap.min_identity`, `AmrHeatmap.title`, `AmrHeatmap.collection` (préfixe
des fichiers de sortie).

## Passage à Terra

Mêmes clés, mais les chemins de fichiers deviennent des `gs://…` (voir
`inputs.example.json` à la racine). Le `docker` doit alors être **poussé** sur
Docker Hub (public).
