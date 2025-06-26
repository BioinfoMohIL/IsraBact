version 1.0

task mst_graph {
    input {
        File alleles_tsv
        String output_png
    }

    command <<<
        echo "Generating Minimal Spanning Tree graph..."

        python3 - <<CODE
            import numpy as np
            import pandas as pd
            import networkx as nx
            import matplotlib.pyplot as plt
            from scipy.spatial.distance import pdist, squareform
            from scipy.sparse.csgraph import minimum_spanning_tree

            alleles = pd.read_csv("~{alleles_tsv}", sep="\\t", header=0, index_col=0)
            alleles.replace("", np.nan, inplace=True)
            alleles = alleles.apply(pd.to_numeric, errors='coerce')

            def hamming_ignore_nan(u, v):
                mask = ~np.logical_or(np.isnan(u), np.isnan(v))
                return np.sum(u[mask] != v[mask]) if np.sum(mask) > 0 else np.nan

            dist_array = pdist(alleles.values, metric=hamming_ignore_nan)
            dist_matrix = squareform(dist_array)
            dist_matrix = np.nan_to_num(dist_matrix, nan=0)

            mst_sparse = minimum_spanning_tree(dist_matrix)
            mst_dense = mst_sparse.toarray().astype(float)
            G = nx.from_numpy_array(mst_dense)
            G = nx.relabel_nodes(G, {i: name for i, name in enumerate(alleles.index)})

            pos = nx.kamada_kawai_layout(G, weight='weight')

            min_length = 0.15
            for u, v in G.edges():
                dx = pos[v][0] - pos[u][0]
                dy = pos[v][1] - pos[u][1]
                dist = np.hypot(dx, dy)
                if dist < min_length and dist > 0:
                    scale = min_length / dist
                    mx, my = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
                    pos[u] = (mx + (pos[u][0] - mx) * scale, my + (pos[u][1] - my) * scale)
                    pos[v] = (mx + (pos[v][0] - mx) * scale, my + (pos[v][1] - my) * scale)

            edges = G.edges(data=True)
            weights = [d['weight'] for (u, v, d) in edges]
            edge_widths = [max(1.0, 10.0 / (w + 1)) for w in weights]

            plt.figure(figsize=(14, 14))
            nx.draw_networkx_nodes(G, pos, node_size=1000, node_color="#ADD8E6", linewidths=0)
            nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
            nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color='gray')
            nx.draw_networkx_edge_labels(
                G, pos,
                edge_labels={(u, v): f"{d['weight']:.0f}" for u, v, d in edges},
                font_size=8, font_color='blue'
            )

            plt.title("cgMLST Allelic Distances", fontsize=16)
            plt.axis('off')
            plt.tight_layout()
            plt.savefig("~{output_png}", dpi=300)
            CODE
    >>>

    output {
        File graph = output_png
    }

    runtime {
        docker: "bioinfomoh/data_analysis_pytools:1"
    }
}
