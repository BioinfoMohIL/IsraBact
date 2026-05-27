"""
draw_mst_network.py
───────────────────
Génère un graphe MST style GrapeTree (nœuds circulaires colorés, distances
sur les branches) depuis un fichier .nwk et un fichier de métadonnées optionnel.

Usage :
    # Sans métadonnées (nœuds gris uniformes)
    python draw_mst_network.py --input tree.nwk --output mst.png

    # Avec métadonnées (couleur par année ou autre colonne)
    python draw_mst_network.py --input tree.nwk --meta metadata.tsv \
                               --color-col year --output mst.png

Format metadata.tsv (TSV, séparateur tabulation) :
    sample_id   year    country   ...
    CA110       2019    Canada
    CA114       2022    Canada

Dépendances :
    pip install biopython matplotlib networkx numpy
"""

import argparse
import csv
import io
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import networkx as nx
import numpy as np
from Bio import Phylo


# ── Palette années (calquée sur GrapeTree) ────────────────────────────────────
YEAR_PALETTE = {
    '2018': '#FDAE61',  # orange clair
    '2019': '#4575B4',  # bleu foncé
    '2020': '#A8D96F',  # vert clair
    '2021': '#D73027',  # rouge
    '2022': '#F46D43',  # orange
    '2023': '#1A9850',  # vert foncé
    '2024': '#ABD9E9',  # bleu clair
}
DEFAULT_NODE_COLOR = '#AAAAAA'
EDGE_COLOR         = '#222222'
BG_COLOR           = '#FFFFFF'
LABEL_COLOR        = '#FFFFFF'
DIST_COLOR         = '#888888'


# ── Helpers ───────────────────────────────────────────────────────────────────
def newick_to_graph(tree) -> nx.Graph:
    G = nx.Graph()

    def walk(clade, parent_name=None):
        name = clade.name if clade.name else f"_internal_{id(clade)}"
        bl   = clade.branch_length or 0
        G.add_node(name, is_leaf=len(clade.clades) == 0)
        if parent_name is not None:
            G.add_edge(parent_name, name, weight=bl)
        for child in clade.clades:
            walk(child, name)

    walk(tree.root)
    return G

def contract_internal_nodes(G: nx.Graph) -> nx.Graph:
    """Supprime les nœuds internes en fusionnant leurs edges."""
    G = G.copy()
    internal = [n for n, d in G.nodes(data=True) 
                if n.startswith('_internal_') or not d.get('is_leaf', True)]
    
    for node in internal:
        neighbors = list(G.neighbors(node))
        if len(neighbors) == 2:
            # Nœud de passage → fusionner les deux edges
            u, v = neighbors
            w = G[node][u]['weight'] + G[node][v]['weight']
            G.add_edge(u, v, weight=w)
        G.remove_node(node)
    
    return G
def load_metadata(path: str, id_col: str = 'sample_id') -> dict:
    """Charge un TSV de métadonnées. Retourne {sample_id: {col: val}}."""
    meta = {}
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = row.get(id_col, '').strip()
            if key:
                meta[key] = {k: v.strip() for k, v in row.items()}
    return meta


def get_color_mapping(G: nx.Graph, meta: dict, color_col: str):
    """
    Retourne (node_colors, palette) où node_colors = {node: hex_color}.
    Gère les valeurs discrètes (année, pays…).
    """
    # Collecter toutes les valeurs uniques
    values = sorted(set(
        meta[n][color_col]
        for n in G.nodes
        if n in meta and color_col in meta[n] and meta[n][color_col]
    ))

    if not values:
        return {n: DEFAULT_NODE_COLOR for n in G.nodes}, {}

    # Utiliser la palette prédéfinie si dispo (années), sinon générer
    palette = {}
    for v in values:
        if v in YEAR_PALETTE:
            palette[v] = YEAR_PALETTE[v]

    # Couleurs manquantes → colormap tab20
    missing = [v for v in values if v not in palette]
    if missing:
        cmap   = plt.cm.get_cmap('tab20', max(len(missing), 1))
        taken  = set(palette.values())
        ci     = 0
        for v in missing:
            while mcolors.to_hex(cmap(ci)) in taken:
                ci += 1
            palette[v] = mcolors.to_hex(cmap(ci))
            ci += 1

    node_colors = {}
    for n in G.nodes:
        if n in meta and color_col in meta[n] and meta[n][color_col]:
            node_colors[n] = palette.get(meta[n][color_col], DEFAULT_NODE_COLOR)
        else:
            node_colors[n] = DEFAULT_NODE_COLOR

    return node_colors, palette


def node_size(G: nx.Graph, n: str, meta: dict, count_col: str | None) -> float:
    """
    Taille du nœud proportionnelle au nombre de samples si count_col fourni,
    sinon proportionnelle au degré (comme GrapeTree).
    """
    base = 1200
    if count_col and n in meta:
        try:
            count = int(meta[n].get(count_col, 1))
            return base * math.sqrt(count)
        except ValueError:
            pass
    deg = G.degree(n)
    return base * max(0.6, math.sqrt(deg))


# ── Draw ──────────────────────────────────────────────────────────────────────

def draw_mst(
    G: nx.Graph,
    node_colors: dict,
    palette: dict,
    color_col: str | None,
    meta: dict,
    output_path: str,
    title: str,
    dpi: int,
    count_col: str | None = None,
):
    fig, ax = plt.subplots(figsize=(14, 10))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)
    ax.axis('off')

    # Layout — spring avec poids inversés (branches courtes = proches)
    # On utilise les distances comme distance cible
    weights = nx.get_edge_attributes(G, 'weight')
    max_w   = max(weights.values(), default=1) or 1

    # Convertir poids en longueur pour spring layout
    for u, v in G.edges():
        G[u][v]['spring_len'] = G[u][v]['weight'] / max_w

    pos = nx.spring_layout(
        G,
        weight     = 'spring_len',
        seed       = 42,
        k          = 2.5 / math.sqrt(max(G.number_of_nodes(), 1)),
        iterations = 200,
    )

    # ── Edges ─────────────────────────────────────────────────────
    # Épaisseur proportionnelle à 1/distance (branches courtes = plus épaisses)
    for u, v, data in G.edges(data=True):
        w   = data.get('weight', 1)
        lw  = max(0.5, 2.5 - (w / max_w) * 2)
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)],
            edge_color=EDGE_COLOR, width=lw, ax=ax, alpha=0.85
        )

    # ── Distance labels sur les branches ──────────────────────────
    edge_labels = {
        (u, v): str(int(round(d['weight'])))
        for u, v, d in G.edges(data=True)
        if d.get('weight', 0) > 0
    }
    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels     = edge_labels,
        font_size       = 7,
        font_color      = DIST_COLOR,
        font_family     = 'sans-serif',
        bbox            = dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7),
        ax              = ax,
    )

    # ── Nodes ──────────────────────────────────────────────────────
    sizes  = [node_size(G, n, meta, count_col) for n in G.nodes()]
    colors = [node_colors.get(n, DEFAULT_NODE_COLOR) for n in G.nodes()]

    nx.draw_networkx_nodes(
        G, pos,
        node_size  = sizes,
        node_color = colors,
        linewidths = 1.2,
        edgecolors = '#444444',
        ax         = ax,
    )

    # ── Node labels ────────────────────────────────────────────────
    # Afficher seulement les noms de feuilles (pas les nœuds internes node_xxx)
    display_labels = {
        n: n for n in G.nodes()
        if not n.startswith('node_')
    }
    nx.draw_networkx_labels(
        G, pos,
        labels    = display_labels,
        font_size = 7,
        font_color = LABEL_COLOR,
        font_weight = 'bold',
        ax        = ax,
    )

    # ── Légende ────────────────────────────────────────────────────
    if palette and color_col:
        # Compter les samples par valeur si méta disponible
        count_by_val = {}
        for n, m in meta.items():
            v = m.get(color_col, '')
            if v:
                count_by_val[v] = count_by_val.get(v, 0) + 1

        handles = []
        for val in sorted(palette):
            count = count_by_val.get(val, '')
            label = f"{val} [{count}]" if count else val
            handles.append(mpatches.Patch(
                facecolor = palette[val],
                edgecolor = '#444',
                linewidth = 0.8,
                label     = label,
            ))

        legend = ax.legend(
            handles        = handles,
            title          = color_col.capitalize(),
            title_fontsize = 9,
            fontsize       = 8,
            loc            = 'upper right',
            framealpha     = 0.9,
            edgecolor      = '#ddd',
            facecolor      = 'white',
            borderpad      = 0.8,
        )
        legend.get_title().set_fontweight('bold')

    # ── Scale bar (200 allelic differences) ───────────────────────
    # Calculer la longueur en unités de position pour 200 allèles
    scale_alleles = 200
    scale_frac    = scale_alleles / max_w * 0.08  # approximation visuelle
    ax.annotate(
        '', xy=(0.88, 0.03), xytext=(0.88 + scale_frac, 0.03),
        xycoords='axes fraction',
        arrowprops=dict(arrowstyle='<->', color='black', lw=1.2)
    )
    ax.text(0.88 + scale_frac / 2, 0.01, str(scale_alleles),
            transform=ax.transAxes, ha='center', va='top',
            fontsize=8, color='black')

    ax.set_title(title, fontsize=13, fontweight='bold',
                 color='#1B4332', pad=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    print(f"[✓] MST network saved → {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Draw MST network (GrapeTree style) from .nwk + optional metadata'
    )
    parser.add_argument('--input',     '-i', required=False,
                        help='Path to .nwk file (Newick)')
    parser.add_argument('--meta',      '-m', default=None,
                        help='Path to metadata TSV (optional)')
    parser.add_argument('--id-col',          default='sample_id',
                        help='Column name for sample IDs in metadata [default: sample_id]')
    parser.add_argument('--color-col', '-c', default='year',
                        help='Metadata column to color nodes by [default: year]')
    parser.add_argument('--count-col',       default=None,
                        help='Metadata column with sample count (for node size)')
    parser.add_argument('--output',    '-o', default='mst_network.png',
                        help='Output PNG path [default: mst_network.png]')
    parser.add_argument('--title',     '-t', default='Minimum Spanning Tree — cgMLST',
                        help='Plot title')
    parser.add_argument('--dpi',             type=int, default=200,
                        help='Resolution [default: 200]')
    args = parser.parse_args()

    # ── Charger le NWK ────────────────────────────────────────────
    TEST_NWK = (
        "((((126:84,219:11,223:11,217:213,208:103)48_a:48,"
        "(46:174,97:15,118:89,84:26,120:26,211:36,222:36)48_b:43)221_a:48,"
        "(220:0,214:0)221_b:0,(110:1228,176:1209)221_c:596)root_a:538,"
        "((89:58,91:66,111:105,114:158,215:62,175:30,169:115,"
        "(166:42,103:42)169_a:377,117:115)174_a:140,"
        "(63:66)174_b:140)174_c:538);"
    )

    if args.input:
        with open(args.input) as f:
            nwk = f.read().strip()
    else:
        print("[info] No --input provided, using built-in test tree.")
        nwk = TEST_NWK

    tree = Phylo.read(io.StringIO(nwk), 'newick')
    G    = newick_to_graph(tree)
    G = contract_internal_nodes(G)
    print(f"[info] Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # ── Charger les métadonnées (optionnel) ───────────────────────
    meta = {}
    if args.meta:
        meta = load_metadata(args.meta, id_col=args.id_col)
        print(f"[info] Metadata loaded: {len(meta)} samples")

    # ── Générer des métadonnées fictives si absentes + test ───────
    if not meta:
        years = ['2019', '2020', '2021', '2022', '2023', '2024']
        rng   = np.random.default_rng(42)
        for n in G.nodes():
            if not n.startswith('node_'):
                meta[n] = {'year': rng.choice(years), 'sample_id': n}
        print("[info] No metadata provided — using random years for demo.")

    # ── Couleurs ──────────────────────────────────────────────────
    node_colors, palette = get_color_mapping(G, meta, args.color_col)

    # ── Draw ──────────────────────────────────────────────────────
    draw_mst(
        G            = G,
        node_colors  = node_colors,
        palette      = palette,
        color_col    = args.color_col,
        meta         = meta,
        output_path  = args.output,
        title        = args.title,
        dpi          = args.dpi,
        count_col    = args.count_col,
    )


if __name__ == '__main__':
    main()