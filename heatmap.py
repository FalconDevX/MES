import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np

def plot_heatmap(grid, T):
    T = T.flatten()
    patches = []
    colors = []

    nodes_map = {n.id: n for n in grid.nodes}

    for element in grid.elements:
        coords = []
        temps = []

        for nid in element.nodes_id:
            node = nodes_map[nid]
            coords.append([node.x, node.y])
            temps.append(T[nid - 1]) 

        patches.append(Polygon(coords, closed=True))
        colors.append(np.mean(temps)) 

    fig, ax = plt.subplots(figsize=(6, 5))
    pc = PatchCollection(
        patches,
        cmap="inferno",
        edgecolor="black",
        linewidths=0.1
    )
    pc.set_array(np.array(colors))
    ax.add_collection(pc)
    ax.autoscale()

    plt.colorbar(pc, ax=ax, label="Temperature [°C]")
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Heatmap")

    plt.tight_layout()
    plt.show()