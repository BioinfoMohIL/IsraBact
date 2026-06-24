#!/usr/bin/env python3
"""
plot_heatmap.py - "Fig 3e" style figure: phylogeny on the left, a
resistome/virulome heatmap on the right, optionally annotated by ST/species.

Inputs:
  --matrix       TSV matrix (samples x genes) produced by parse_amr.py
  --genes        column annotation (gene, category, drug_class)    [optional]
  --tree         Newick tree (.nwk); sets the row order            [optional]
  --metadata     TSV with a 'sample' column + columns to annotate  [optional]
  --annotate     metadata columns to display (e.g. ST species)     [optional]

Without --tree: samples are sorted alphabetically (no dendrogram, to keep
dependencies minimal and the Docker image light).

Dependencies: pandas, numpy, matplotlib, biopython. No Qt, scipy or R.
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

try:
    from Bio import Phylo
    HAVE_BIO = True
except Exception:
    HAVE_BIO = False

CAT_COLORS = {
    "resistance": "#c0392b",
    "virulence": "#8e44ad",
    "stress": "#d68910",
    "plasmid": "#16a085",
    "other": "#7f8c8d",
}


# --- tree: positions computed manually to line up exactly with the heatmap ---
def layout_tree(tree):
    """Return x (cumulative depth) and y (tip = index 0..N-1 from top) per clade."""
    tree.ladderize()
    tips = tree.get_terminals()
    y = {}
    for i, t in enumerate(tips):
        y[t] = i  # 0 en haut -> N-1 en bas (origin 'upper' cote heatmap)

    # profondeur (somme des longueurs de branche depuis la racine)
    x = {}
    def set_depth(clade, depth):
        bl = clade.branch_length if clade.branch_length else 0.0
        d = depth + bl
        x[clade] = d
        for c in clade.clades:
            set_depth(c, d)
    set_depth(tree.root, 0.0)

    # y of internal nodes = mean of children
    def set_y(clade):
        if clade.is_terminal():
            return y[clade]
        ys = [set_y(c) for c in clade.clades]
        y[clade] = sum(ys) / len(ys)
        return y[clade]
    set_y(tree.root)

    order = [t.name for t in tips]
    return x, y, order


def draw_tree(ax, tree):
    x, y, order = layout_tree(tree)
    for clade in tree.find_clades():
        if not clade.clades:
            continue
        # vertical bar linking the children
        child_ys = [y[c] for c in clade.clades]
        ax.plot([x[clade], x[clade]], [min(child_ys), max(child_ys)],
                color="#2c3e50", lw=1.0)
        # horizontal branches to each child
        for c in clade.clades:
            ax.plot([x[clade], x[c]], [y[c], y[c]], color="#2c3e50", lw=1.0)
    ax.set_xlim(-0.02 * max(x.values() or [1]), max(x.values() or [1]) * 1.05)
    ax.set_ylim(len(order) - 0.5, -0.5)  # haut = index 0
    ax.axis("off")
    return order


def main():
    ap = argparse.ArgumentParser(description="Phylogenie + heatmap AMR (style Fig 3e)")
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--genes", default=None)
    ap.add_argument("--tree", default=None)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--annotate", nargs="*", default=[])
    ap.add_argument("--title", default="Resistome / virulome")
    ap.add_argument("--out", default="amr_heatmap.pdf")
    ap.add_argument("--png", default=None, help="additional PNG path (optional)")
    ap.add_argument("--svg", default=None, help="additional SVG path (optional)")
    ap.add_argument("--html", default=None,
                    help="self-contained HTML (inline SVG + PNG/PDF/SVG download buttons)")
    args = ap.parse_args()

    mat = pd.read_csv(args.matrix, sep="\t", index_col=0)

    # row order
    order = None
    tree = None
    if args.tree:
        if not HAVE_BIO:
            sys.exit("biopython requis pour --tree")
        tree = Phylo.read(args.tree, "newick")
        _, _, tip_order = layout_tree(tree)
        present = [s for s in tip_order if s in mat.index]
        missing = [s for s in tip_order if s not in mat.index]
        if missing:
            sys.stderr.write(f"[warn] {len(missing)} feuilles de l'arbre absentes de la matrice: "
                             f"{missing[:5]}{'...' if len(missing) > 5 else ''}\n")
        extra = [s for s in mat.index if s not in tip_order]
        if extra:
            sys.stderr.write(f"[warn] {len(extra)} echantillons hors de l'arbre (ignores pour l'ordre)\n")
        order = present
    else:
        order = sorted(mat.index)

    mat = mat.reindex(order)
    n = len(order)
    if n == 0:
        sys.exit("Aucun echantillon a tracer (intersection arbre/matrice vide?).")

    # gene annotation -> order columns by category then name
    if args.genes:
        gann = pd.read_csv(args.genes, sep="\t")
        cat_of = dict(zip(gann["gene"], gann["category"]))
    else:
        cat_of = {}
    cat_of = {g: cat_of.get(g, "resistance") for g in mat.columns}
    col_order = sorted(mat.columns, key=lambda g: (cat_of[g], g))
    mat = mat[col_order]
    col_cats = [cat_of[g] for g in col_order]

    # metadata / side annotations
    annot_cols = []
    if args.metadata and args.annotate:
        meta = pd.read_csv(args.metadata, sep="\t", dtype=str).set_index("sample")
        for col in args.annotate:
            if col in meta.columns:
                annot_cols.append((col, meta[col].reindex(order).fillna("NA")))

    # --- mise en page ---
    has_tree = tree is not None
    ncols = mat.shape[1]
    width_ratios = ([1.6] if has_tree else []) + [0.12] * len(annot_cols) + [max(4, ncols * 0.16)]
    fig_w = 2.2 + 0.16 * ncols + 0.6 * len(annot_cols)
    fig_h = max(3.0, 0.28 * n + 1.4)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(1, len(width_ratios), width_ratios=width_ratios, wspace=0.04)

    ci = 0
    if has_tree:
        ax_tree = fig.add_subplot(gs[0, ci]); ci += 1
        draw_tree(ax_tree, tree)
        ax_tree.set_title("Phylogenie", fontsize=9)

    # colonnes d'annotation (ST, espece...)
    for label, series in annot_cols:
        axa = fig.add_subplot(gs[0, ci]); ci += 1
        cats = pd.Categorical(series.values)
        codes = cats.codes.reshape(-1, 1)
        cmap = plt.get_cmap("tab20", max(1, len(cats.categories)))
        axa.imshow(codes, aspect="auto", cmap=cmap, vmin=0, vmax=max(1, len(cats.categories) - 1))
        axa.set_xticks([]); axa.set_yticks([])
        axa.set_xlabel(label, fontsize=8, rotation=90, labelpad=6)

    # main heatmap
    ax = fig.add_subplot(gs[0, ci])
    binary = set(np.unique(mat.values)) <= {0.0, 1.0}
    if binary:
        cmap = ListedColormap(["#ecf0f1", "#2c3e50"])
        ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    else:
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
        cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("% identite", fontsize=8)

    ax.set_yticks(range(n))
    ax.set_yticklabels(order, fontsize=7)
    if not has_tree and not annot_cols:
        ax.yaxis.tick_left()
    else:
        ax.yaxis.set_visible(False)
    ax.set_xticks(range(ncols))
    ax.set_xticklabels(col_order, fontsize=6, rotation=90)
    ax.set_title(args.title, fontsize=10)

    # gene-category band above the heatmap
    used_cats = []
    for j, c in enumerate(col_cats):
        ax.add_patch(plt.Rectangle((j - 0.5, -1.2), 1, 0.8, clip_on=False,
                                   color=CAT_COLORS.get(c, CAT_COLORS["other"])))
        if c not in used_cats:
            used_cats.append(c)
    handles = [Patch(color=CAT_COLORS.get(c, CAT_COLORS["other"]), label=c) for c in used_cats]
    if binary:
        handles += [Patch(color="#2c3e50", label="present"),
                    Patch(color="#ecf0f1", label="absent")]
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=7, frameon=False)

    fig.suptitle("AMR resistome/virulome heatmap", fontsize=12, y=0.995)
    fig.savefig(args.out, bbox_inches="tight", dpi=200)
    if args.png:
        fig.savefig(args.png, bbox_inches="tight", dpi=200)
    if args.svg:
        fig.savefig(args.svg, bbox_inches="tight")
    if args.html:
        write_html(fig, args.html, args.title)
    extra = [p for p in (args.png, args.svg, args.html) if p]
    sys.stderr.write(f"[ok] figure -> {args.out}"
                     + (" + " + " + ".join(extra) if extra else "") + "\n")


def write_html(fig, path, title):
    """Self-contained HTML: crisp inline SVG + buttons to download PNG/PDF/SVG."""
    import base64, io
    svg_buf = io.StringIO()
    fig.savefig(svg_buf, format="svg", bbox_inches="tight")
    svg = svg_buf.getvalue()
    svg = svg[svg.find("<svg"):]  # strip XML/doctype preamble

    def b64(fmt):
        b = io.BytesIO()
        fig.savefig(b, format=fmt, bbox_inches="tight",
                    dpi=(300 if fmt == "png" else None))
        return base64.b64encode(b.getvalue()).decode()

    png64, pdf64, svg64 = b64("png"), b64("pdf"), base64.b64encode(svg.encode()).decode()
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title} - AMR heatmap</title>
<style>
  :root {{ --ink:#1f2933; --bg:#ffffff; --accent:#2c3e50; }}
  body {{ margin:0; font-family:-apple-system,Segoe UI,Roboto,sans-serif;
          color:var(--ink); background:#eef1f4; }}
  header {{ padding:18px 24px; background:var(--bg); border-bottom:1px solid #e2e6ea;
            display:flex; align-items:center; gap:16px; flex-wrap:wrap; }}
  h1 {{ font-size:18px; margin:0; font-weight:600; }}
  .btns {{ margin-left:auto; display:flex; gap:8px; }}
  button {{ border:1px solid var(--accent); background:var(--accent); color:#fff;
            padding:8px 14px; border-radius:8px; font-size:13px; cursor:pointer; }}
  button.alt {{ background:#fff; color:var(--accent); }}
  button:hover {{ opacity:.9; }}
  main {{ padding:24px; display:flex; justify-content:center; }}
  .fig {{ background:var(--bg); padding:16px; border-radius:12px;
          box-shadow:0 1px 4px rgba(0,0,0,.08); max-width:100%; overflow:auto; }}
  .fig svg {{ max-width:100%; height:auto; display:block; }}
  footer {{ text-align:center; color:#7b8794; font-size:12px; padding:0 0 24px; }}
</style></head>
<body>
<header>
  <h1>{title}</h1>
  <div class="btns">
    <button onclick="dl('png','image/png','heatmap.png')">Download PNG</button>
    <button class="alt" onclick="dl('pdf','application/pdf','heatmap.pdf')">PDF</button>
    <button class="alt" onclick="dl('svg','image/svg+xml','heatmap.svg')">SVG</button>
  </div>
</header>
<main><div class="fig">{svg}</div></main>
<footer>Generated by plot_heatmap.py - resistome/virulome heatmap</footer>
<script>
  const data = {{ png:"{png64}", pdf:"{pdf64}", svg:"{svg64}" }};
  function dl(fmt, mime, name) {{
    const a = document.createElement('a');
    a.href = 'data:' + mime + ';base64,' + data[fmt];
    a.download = name; document.body.appendChild(a); a.click(); a.remove();
  }}
</script>
</body></html>"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)


if __name__ == "__main__":
    main()