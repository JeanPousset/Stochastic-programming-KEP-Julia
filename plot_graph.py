# %% Finction
import networkx as nx
import matplotlib.pyplot as plt

def read_wmd_edges(filename):
    """
    Lit un fichier .wmd et retourne une liste d'arêtes (u, v, weight)
    en ignorant les lignes de métadonnées commençant par '#'.
    """
    edges = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # ignore les métadonnées
            parts = line.split(",")
            if len(parts) >= 2:
                u, v = int(parts[0]), int(parts[1])
                w = float(parts[2]) if len(parts) > 2 else 1.0
                edges.append((u, v, w))
    return edges

def plot_graph_from_file(filename, output="graph.png"):
    # Lecture des arêtes
    edges = read_wmd_edges(filename)
    
    # Création du graphe orienté
    G = nx.DiGraph()
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)
    
    # Mise en page du graphe
    pos = nx.spring_layout(G, seed=42, k=0.8/len(G.nodes())**0.5)
    
    # Tracé
    plt.figure(figsize=(8, 6))
    nx.draw_networkx_nodes(G, pos, node_size=500, node_color="#99ccff", edgecolors="black")
    nx.draw_networkx_edges(G, pos, arrowstyle="->", arrowsize=15, width=1)
    nx.draw_networkx_labels(G, pos, font_size=9, font_color="black")
    
    plt.title(f"Graphe issu de {filename}", fontsize=12)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    print(f"✅ Graphe sauvegardé dans {output}")


# %% 
plot_graph_from_file("../KEP_001.wmd", "G_001.png")
