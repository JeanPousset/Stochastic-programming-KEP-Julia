# %% Finction
import networkx as nx
import matplotlib.pyplot as plt

def read_wmd_edges(filename):
    edges = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) >= 2:
                u, v = int(parts[0]), int(parts[1])
                w = float(parts[2]) if len(parts) > 2 else 1.0
                edges.append((u, v, w))
    return edges

def plot_graph_from_file(filename, output="graph.png"):
    edges = read_wmd_edges(filename)
    G = nx.DiGraph()
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)

    # Choix du layout (testé pour graphes denses)
    pos = nx.kamada_kawai_layout(G)

    # Taille automatique
    n = len(G.nodes())
    plt.figure(figsize=(max(10, n/1.5), max(8, n/1.5)))

    nx.draw_networkx_nodes(G, pos, node_size=600, node_color="#a7c7e7", edgecolors="black")
    nx.draw_networkx_edges(G, pos, arrowstyle="->", arrowsize=15, width=1)
    nx.draw_networkx_labels(G, pos, font_size=10)

    plt.title(f"Graphe issu de {filename}", fontsize=12)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    print(f"✅ Graphe sauvegardé dans {output}")


# %% 
plot_graph_from_file("../KEP_001.wmd", "G_001.png")


# %%
