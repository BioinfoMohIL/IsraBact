#!/usr/bin/env python3
"""
plot_heatmap.py - AMRViz "Fig 3e" style figure: a rectangular cladogram on the
left aligned to a presence/absence heatmap, where each present cell is colored
by gene category (resistance = blue, virulence = red, ...) and absent = white.

Inputs:
  --matrix       TSV matrix (samples x genes) from parse_amr.py
  --genes        column annotation (gene, category, drug_class)   [optional]
  --tree         Newick tree (.nwk); sets the row order           [optional]
  --metadata     TSV with a 'sample' column + columns to annotate [optional]
  --annotate     metadata columns to display                      [optional]

Outputs: PDF, and optionally PNG / SVG / a self-contained HTML.
No titles are drawn on the figure.
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

try:
    from Bio import Phylo
    HAVE_BIO = True
except Exception:
    HAVE_BIO = False

# present-cell color per category; absent cells are white
CELL = {
    "resistance": "#2c5f8a",   # blue
    "virulence":  "#a93226",   # red
    "stress":     "#ca6f1e",   # orange
    "plasmid":    "#148f77",   # teal
    "other":      "#7f8c8d",
}
BAND_LABEL = {
    "resistance": "AMR", "virulence": "Virulence",
    "stress": "Stress", "plasmid": "Plasmid", "other": "Other",
}
# left-to-right block order on the heatmap
CAT_ORDER = ["virulence", "resistance", "stress", "plasmid", "other"]


# --- rectangular cladogram: tips aligned right, positions match heatmap rows ---
def layout_cladogram(tree):
    tree.ladderize()
    tips = tree.get_terminals()
    y = {t: i for i, t in enumerate(tips)}          # 0..N-1 from top

    # height = number of edges on the longest path down to a leaf
    height = {}
    def get_height(c):
        height[c] = 0 if c.is_terminal() else 1 + max(get_height(ch) for ch in c.clades)
        return height[c]
    root_h = get_height(tree.root)

    x = {c: root_h - height[c] for c in tree.find_clades()}  # tips at x=root_h
    def set_y(c):
        if c.is_terminal():
            return y[c]
        ys = [set_y(ch) for ch in c.clades]
        y[c] = sum(ys) / len(ys)
        return y[c]
    set_y(tree.root)
    return x, y, [t.name for t in tips], root_h


def draw_cladogram(ax, tree):
    x, y, order, root_h = layout_cladogram(tree)
    for clade in tree.find_clades():
        if not clade.clades:
            continue
        child_ys = [y[c] for c in clade.clades]
        ax.plot([x[clade], x[clade]], [min(child_ys), max(child_ys)],
                color="#2c3e50", lw=1.0)
        for c in clade.clades:
            ax.plot([x[clade], x[c]], [y[c], y[c]], color="#2c3e50", lw=1.0)
    ax.set_xlim(-0.2, root_h + 0.2)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.axis("off")
    return order


def cladogram_segments(tree):
    """Return tree line segments in (x in [0,root_h], y in [0,n-1]) space for JS."""
    x, y, order, root_h = layout_cladogram(tree)
    segs = []
    for clade in tree.find_clades():
        if not clade.clades:
            continue
        cy = [y[c] for c in clade.clades]
        segs.append([x[clade], min(cy), x[clade], max(cy)])
        for c in clade.clades:
            segs.append([x[clade], y[c], x[c], y[c]])
    return segs, root_h, order



def category_blocks(col_cats):
    """Contiguous (category, start, end) column ranges."""
    blocks, j = [], 0
    while j < len(col_cats):
        k = j
        while k + 1 < len(col_cats) and col_cats[k + 1] == col_cats[j]:
            k += 1
        blocks.append((col_cats[j], j, k))
        j = k + 1
    return blocks


def draw_cells(ax, present, col_cats, border="#eef1f3", lw=0.4):
    """Vector cells: present -> category color, absent -> white. Thin borders."""
    n, m = present.shape
    for j in range(m):
        color = CELL.get(col_cats[j], CELL["other"])
        for i in range(n):
            ax.add_patch(plt.Rectangle(
                (j - 0.5, i - 0.5), 1, 1,
                facecolor=(color if present[i, j] else "#ffffff"),
                edgecolor=border, linewidth=lw))


def draw_band(ax, blocks, n, gap, bandh):
    """Colored category band above the data (with a gap) + block separators."""
    ytop = -0.5 - gap - bandh
    for cat, start, end in blocks:
        ax.add_patch(plt.Rectangle((start - 0.5, ytop), end - start + 1, bandh,
                     facecolor=CELL.get(cat, CELL["other"]), edgecolor="none",
                     clip_on=False))
        ax.text((start + end) / 2, ytop + bandh / 2, BAND_LABEL.get(cat, cat),
                ha="center", va="center", color="white", fontsize=8,
                fontweight="bold", clip_on=False)
    for cat, start, end in blocks[:-1]:
        ax.plot([end + 0.5, end + 0.5], [-0.5, n - 0.5], color="white", lw=1.4)


def draw_names(ax, order, ylim):
    """Sample names in their own column (no overlap with the tree)."""
    ax.set_xlim(0, 1); ax.set_ylim(*ylim); ax.axis("off")
    for i, s in enumerate(order):
        ax.text(0.96, i, str(s), ha="right", va="center",
                fontsize=7, color="#3e4c59")


def main():
    ap = argparse.ArgumentParser(description="Cladogram + AMR/virulence heatmap (Fig 3e style)")
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--genes", default=None)
    ap.add_argument("--tree", default=None)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--annotate", nargs="*", default=[])
    ap.add_argument("--title", default="Resistome / virulome",
                    help="used only as the HTML page title, never drawn on the figure")
    ap.add_argument("--out", default="amr_heatmap.pdf")
    ap.add_argument("--png", default=None)
    ap.add_argument("--svg", default=None)
    ap.add_argument("--html", default=None)
    args = ap.parse_args()

    mat = pd.read_csv(args.matrix, sep="\t", index_col=0)

    # row order from the tree (cladogram tip order), else alphabetical
    tree = None
    if args.tree:
        if not HAVE_BIO:
            sys.exit("biopython required for --tree")
        tree = Phylo.read(args.tree, "newick")
        _, _, tip_order, _ = layout_cladogram(tree)
        order = [s for s in tip_order if s in mat.index]
        missing = [s for s in tip_order if s not in mat.index]
        if missing:
            sys.stderr.write(f"[warn] {len(missing)} tree leaves absent from matrix: {missing[:5]}\n")
    else:
        order = sorted(mat.index)
    mat = mat.reindex(order)
    n = len(order)
    if n == 0:
        sys.exit("No samples to plot (empty tree/matrix intersection?).")

    # category per gene, then order columns by category block + name
    cat_of = {}
    if args.genes:
        g = pd.read_csv(args.genes, sep="\t")
        cat_of = dict(zip(g["gene"], g["category"]))
    cat_of = {gene: cat_of.get(gene, "resistance") for gene in mat.columns}
    rank = {c: i for i, c in enumerate(CAT_ORDER)}
    col_order = sorted(mat.columns, key=lambda gene: (rank.get(cat_of[gene], 99), gene))
    mat = mat[col_order]
    col_cats = [cat_of[gene] for gene in col_order]
    ncols = len(col_order)

    # optional metadata side columns
    annot_cols = []
    if args.metadata and args.annotate:
        meta = pd.read_csv(args.metadata, sep="\t", dtype=str).set_index("sample")
        for col in args.annotate:
            if col in meta.columns:
                annot_cols.append((col, meta[col].reindex(order).fillna("NA")))

    present = mat.values > 0
    blocks = category_blocks(col_cats)
    has_tree = tree is not None

    # band geometry (in row units) with a clear gap above the data
    bandh, gap = 0.9, 0.5
    ytop = -0.5 - gap - bandh
    ylim = (n - 0.5, ytop - 0.15)

    # --- layout: [tree] [names] [annot...] [heatmap] ---
    name_chars = max((len(str(s)) for s in order), default=4)
    widths = ([1.3] if has_tree else []) + [0.06 + 0.018 * name_chars] \
        + [0.10] * len(annot_cols) + [max(4, ncols * 0.22)]
    fig_w = 1.4 + 0.22 * ncols + 0.5 * len(annot_cols) + (1.5 if has_tree else 0) \
        + 0.09 * name_chars
    fig_h = max(3.0, 0.30 * n + 1.4)
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(1, len(widths), width_ratios=widths, wspace=0.02)

    ci = 0
    if has_tree:
        ax_tree = fig.add_subplot(gs[0, ci]); ci += 1
        draw_cladogram(ax_tree, tree)
        ax_tree.set_ylim(*ylim)

    ax_names = fig.add_subplot(gs[0, ci]); ci += 1
    draw_names(ax_names, order, ylim)

    for label, series in annot_cols:
        axa = fig.add_subplot(gs[0, ci]); ci += 1
        cats = pd.Categorical(series.values)
        cmap = plt.get_cmap("tab20", max(1, len(cats.categories)))
        for i, code in enumerate(cats.codes):
            axa.add_patch(plt.Rectangle((0, i - 0.5), 1, 1,
                          facecolor=cmap(code), edgecolor="white", linewidth=0.4))
        axa.set_xlim(0, 1); axa.set_ylim(*ylim); axa.axis("off")
        axa.text(0.5, ytop + bandh / 2, label, ha="center", va="center",
                 rotation=90, fontsize=7, color="#3e4c59")

    ax = fig.add_subplot(gs[0, ci])
    draw_cells(ax, present, col_cats)
    draw_band(ax, blocks, n, gap, bandh)
    ax.set_xlim(-0.5, ncols - 0.5)
    ax.set_ylim(*ylim)
    ax.set_xticks(range(ncols))
    ax.set_xticklabels(col_order, fontsize=7, rotation=90)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    for p in (args.out,) + tuple(x for x in (args.png, args.svg) if x):
        fig.savefig(p, bbox_inches="tight", pad_inches=0.2, facecolor="white",
                    dpi=300 if str(p).endswith("png") else None)

    if args.html:
        data_blocks = [{"cat": c, "start": s, "end": e,
                        "label": BAND_LABEL.get(c, c)} for c, s, e in blocks]
        segs, root_h, _ = cladogram_segments(tree) if tree is not None else ([], 0, order)
        data = {
            "samples": list(order),
            "genes": list(col_order),
            "cats": list(col_cats),
            "present": (mat.values > 0).astype(int).tolist(),
            "colors": CELL,
            "blocks": data_blocks,
            "tree": {"segments": segs, "root_h": root_h},
        }
        write_html(args.html, args.title, fig, data)

    extra = [p for p in (args.png, args.svg, args.html) if p]
    sys.stderr.write(f"[ok] figure -> {args.out}" + (" + " + " + ".join(extra) if extra else "") + "\n")


def write_html(path, title, fig, data):
    """Self-contained interactive HTML (vanilla JS SVG) + publication downloads."""
    import base64, io, json

    def b64(fmt):
        b = io.BytesIO()
        fig.savefig(b, format=fmt, bbox_inches="tight", dpi=(300 if fmt == "png" else None))
        return base64.b64encode(b.getvalue()).decode()

    assets = {"png": b64("png"), "pdf": b64("pdf"), "svg": b64("svg")}
    head = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>" + str(title) + "</title>"
        "<script>const DATA=" + json.dumps(data) + ";const ASSETS=" + json.dumps(assets)
        + ";const TITLE=" + json.dumps(str(title)) + ";</script>"
    )
    open(path, "w", encoding="utf-8").write(head + HTML_BODY)


HTML_BODY = '\n<style>\n *{box-sizing:border-box}\n body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;background:#f4f6f8;color:#1f2933}\n header{position:sticky;top:0;z-index:5;display:flex;align-items:center;gap:14px;padding:10px 18px;background:#fff;border-bottom:1px solid #e2e6ea}\n h1{font-size:14px;font-weight:600;margin:0;color:#3e4c59}\n .btns{margin-left:auto;display:flex;gap:8px}\n button{border:1px solid #2c3e50;background:#2c3e50;color:#fff;padding:7px 13px;border-radius:7px;font-size:13px;cursor:pointer}\n button.alt{background:#fff;color:#2c3e50} button:hover{opacity:.92}\n main{padding:16px}\n .wrap{background:#fff;border-radius:10px;box-shadow:0 1px 4px rgba(0,0,0,.07);padding:14px;display:inline-block;overflow:auto;max-width:100%}\n svg{display:block;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}\n .tip{position:fixed;pointer-events:none;background:#1f2933;color:#fff;padding:6px 9px;border-radius:6px;font-size:12px;line-height:1.35;opacity:0;transition:opacity .07s;z-index:9;box-shadow:0 2px 8px rgba(0,0,0,.25)}\n .tip b{font-weight:600}\n</style>\n<header>\n  <h1 id="ttl"></h1>\n  <div class="btns">\n    <button onclick="dl(\'png\',\'image/png\',\'heatmap.png\')">Download PNG</button>\n    <button class="alt" onclick="dl(\'svg\',\'image/svg+xml\',\'heatmap.svg\')">SVG</button>\n    <button class="alt" onclick="dl(\'pdf\',\'application/pdf\',\'heatmap.pdf\')">PDF</button>\n  </div>\n</header>\n<main><div class="wrap" id="wrap"></div></main>\n<div class="tip" id="tip"></div>\n<script>\nconst SVGNS="http://www.w3.org/2000/svg";\nconst CELL=22, BAND=24, PADT=8, PADL=10, CHARW=6.4, FONT=11;\ndocument.getElementById(\'ttl\').textContent=TITLE;\nfunction el(t,a){const e=document.createElementNS(SVGNS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}\nfunction maxlen(a){return a.reduce((m,s)=>Math.max(m,(\'\'+s).length),0);}\nfunction R(v){return Math.round(v);}\n\nfunction render(){\n  const S=DATA.samples,G=DATA.genes,C=DATA.cats,P=DATA.present,COL=DATA.colors,B=DATA.blocks;\n  const n=S.length,m=G.length;\n  const nameW=R(maxlen(S)*CHARW+12), geneH=R(maxlen(G)*CHARW+12);\n  const treeW=(DATA.tree&&DATA.tree.segments.length)?150:0;\n  const gx0=PADL+treeW+nameW, gy0=PADT+BAND;\n  const W=gx0+m*CELL+PADL, H=gy0+n*CELL+geneH+PADT;\n  const svg=el(\'svg\',{width:W,height:H,viewBox:`0 0 ${W} ${H}`,xmlns:SVGNS,id:\'heat\',\n                      \'text-rendering\':\'geometricPrecision\'});\n  svg.appendChild(el(\'rect\',{x:0,y:0,width:W,height:H,fill:\'#ffffff\'}));\n\n  // rectangular cladogram: all tips aligned at the same level\n  if(treeW){\n    const rh=DATA.tree.root_h||1;\n    const tx=x=>R(PADL+(rh?(x/rh)*(treeW-18):treeW/2));\n    const ty=y=>R(gy0+y*CELL+CELL/2);\n    const g=el(\'g\',{stroke:\'#2c3e50\',\'stroke-width\':1.1,\'shape-rendering\':\'crispEdges\',fill:\'none\'});\n    for(const s of DATA.tree.segments)\n      g.appendChild(el(\'line\',{x1:tx(s[0]),y1:ty(s[1]),x2:tx(s[2]),y2:ty(s[3])}));\n    svg.appendChild(g);\n  }\n  // sample names\n  for(let i=0;i<n;i++){\n    const t=el(\'text\',{x:gx0-5,y:gy0+i*CELL+CELL/2+3.5,\'text-anchor\':\'end\',\'font-size\':FONT,fill:\'#3e4c59\'});\n    t.textContent=S[i]; svg.appendChild(t);\n  }\n  // category band (one colored bar per gene block) + block separators\n  for(const b of B){\n    const x=gx0+b.start*CELL, w=(b.end-b.start+1)*CELL;\n    svg.appendChild(el(\'rect\',{x:x,y:PADT,width:w,height:BAND-4,fill:COL[b.cat]||\'#7f8c8d\',\n                               \'shape-rendering\':\'crispEdges\'}));\n    const t=el(\'text\',{x:x+w/2,y:PADT+(BAND-4)/2+4,\'text-anchor\':\'middle\',\'font-size\':11,\n                       \'font-weight\':600,fill:\'#fff\'}); t.textContent=b.label; svg.appendChild(t);\n    if(b.end+1<m) svg.appendChild(el(\'line\',{x1:gx0+(b.end+1)*CELL,y1:gy0,\n        x2:gx0+(b.end+1)*CELL,y2:gy0+n*CELL,stroke:\'#fff\',\'stroke-width\':2}));\n  }\n  // cells: present -> category color, absent -> white\n  const cg=el(\'g\',{}); svg.appendChild(cg);\n  for(let i=0;i<n;i++)for(let j=0;j<m;j++){\n    const on=P[i][j];\n    const r=el(\'rect\',{x:gx0+j*CELL,y:gy0+i*CELL,width:CELL,height:CELL,\n        fill:on?(COL[C[j]]||\'#7f8c8d\'):\'#ffffff\',stroke:\'#dfe4e8\',\'stroke-width\':1,\n        \'shape-rendering\':\'crispEdges\'});\n    r.dataset.i=i; r.dataset.j=j; cg.appendChild(r);\n  }\n  // gene names (vertical)\n  for(let j=0;j<m;j++){\n    const cx=gx0+j*CELL+CELL/2, cy=gy0+n*CELL+6;\n    const t=el(\'text\',{x:cx,y:cy,\'font-size\':FONT,fill:\'#3e4c59\',\'text-anchor\':\'end\',\n                       transform:`rotate(-90 ${cx} ${cy})`}); t.textContent=G[j]; svg.appendChild(t);\n  }\n  // hover highlight + tooltip\n  const hl=el(\'rect\',{fill:\'none\',stroke:\'#1f2933\',\'stroke-width\':2,opacity:0}); svg.appendChild(hl);\n  const hr=el(\'rect\',{fill:\'none\',stroke:\'#1f2933\',\'stroke-width\':2,opacity:0}); svg.appendChild(hr);\n  const tip=document.getElementById(\'tip\');\n  cg.addEventListener(\'mousemove\',e=>{\n    if(e.target.tagName!==\'rect\')return;\n    const i=+e.target.dataset.i,j=+e.target.dataset.j;\n    hl.setAttribute(\'x\',gx0+j*CELL);hl.setAttribute(\'y\',gy0);hl.setAttribute(\'width\',CELL);hl.setAttribute(\'height\',n*CELL);hl.setAttribute(\'opacity\',.5);\n    hr.setAttribute(\'x\',gx0);hr.setAttribute(\'y\',gy0+i*CELL);hr.setAttribute(\'width\',m*CELL);hr.setAttribute(\'height\',CELL);hr.setAttribute(\'opacity\',.5);\n    tip.innerHTML=`<b>${S[i]}</b> &middot; ${G[j]}<br>${C[j]} &middot; ${P[i][j]?\'present\':\'absent\'}`;\n    tip.style.left=(e.clientX+14)+\'px\';tip.style.top=(e.clientY+14)+\'px\';tip.style.opacity=1;\n  });\n  cg.addEventListener(\'mouseleave\',()=>{tip.style.opacity=0;hl.setAttribute(\'opacity\',0);hr.setAttribute(\'opacity\',0);});\n  const wrap=document.getElementById(\'wrap\'); wrap.innerHTML=\'\'; wrap.appendChild(svg);\n}\nfunction dl(f,m,n){const a=document.createElement(\'a\');a.href=\'data:\'+m+\';base64,\'+ASSETS[f];a.download=n;document.body.appendChild(a);a.click();a.remove();}\nrender();\n</script></body></html>\n'


if __name__ == "__main__":
    main()