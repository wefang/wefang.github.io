"""Draw the three research-program card images for the about page.

Deliberately schematic: they name each direction's idea (a record read back into a tree,
signals deciding fate over divisions, matching regulation across species) without
reproducing any figure from a proposal or an unpublished manuscript.

    python bin/make_research_art.py      # writes assets/img/research/{lineage,fate,transfer}.png
"""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "research")
os.makedirs(OUT, exist_ok=True)
TEAL, TEAL_L, INK, GREY, LIGHT = "#0F6E6E", "#8CC9C4", "#1B2227", "#9AA7AD", "#E3E8EA"
EDITS = ["#0F6E6E", "#D98E3A", "#6C7FC0", "#C2566B", "#8CB35A"]
W, H = 8.0, 5.0


def canvas():
    fig = plt.figure(figsize=(W, H), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    return fig, ax


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=100, facecolor="white")
    plt.close(fig)


# ---------------- 1. lineage inference: edits accumulate, a dated tree is read back ----------------
fig, ax = canvas()
leaves = np.linspace(0.9, 7.1, 8)
rng = np.random.default_rng(3)


def node(xs, y):
    return (min(xs) + max(xs)) / 2, y


levels = [(0, 4.3), (1, 3.35), (2, 2.45), (3, 1.55)]
groups = [list(range(8))]
coords = {}
for depth, y in levels:
    new = []
    for g in groups:
        x = np.mean([leaves[i] for i in g])
        coords[(depth, tuple(g))] = (x, y)
        if len(g) > 1:
            h = len(g) // 2
            new += [g[:h], g[h:]]
    groups = new
for (depth, g), (x, y) in coords.items():
    if len(g) > 1:
        h = len(g) // 2
        for child in (tuple(g[:h]), tuple(g[h:])):
            cx, cy = coords[(depth + 1, child)] if (depth + 1, child) in coords else (leaves[child[0]], 1.55)
            ax.plot([x, cx, cx], [y, y, cy], color=INK, lw=1.6, solid_capstyle="round")
    ax.add_patch(Circle((x, y), 0.07, color=TEAL, zorder=3))
# barcodes under each leaf: five sites, edits accumulate
for i, x in enumerate(leaves):
    for s in range(5):
        edited = rng.random() < 0.55
        col = EDITS[(i // 2 + s) % 5] if edited else LIGHT
        ax.add_patch(Rectangle((x - 0.3 + s * 0.12, 0.72), 0.11, 0.55, color=col, lw=0))
# time axis
ax.annotate("", xy=(0.42, 1.5), xytext=(0.42, 4.4), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.2))
ax.text(0.14, 2.95, "time", rotation=90, va="center", ha="center", color=GREY, fontsize=20)
save(fig, "lineage.png")

# ---------------- 2. cell fate: a signal gradient decides what dividing cells become ----------------
fig, ax = canvas()
grad = np.linspace(0, 1, 256).reshape(1, -1)
ax.imshow(grad, extent=(0.4, 6.7, 0.35, 0.65), cmap=matplotlib.colors.LinearSegmentedColormap.from_list("g", ["#FFFFFF", TEAL]), aspect="auto")
ax.text(6.8, 0.5, "signal", va="center", ha="left", color=GREY, fontsize=18)
FATE = {0: "#D98E3A", 1: "#6C7FC0", 2: TEAL}


def grow(x, y, span, depth, pos):
    ax.add_patch(Circle((x, y), 0.17, facecolor="white", edgecolor=INK, lw=1.4, zorder=3))
    if depth == 3:
        fate = 0 if pos < 0.33 else (1 if pos < 0.66 else 2)
        ax.add_patch(Circle((x, y), 0.17, facecolor=FATE[fate], edgecolor=INK, lw=1.4, zorder=4))
        return
    for d in (-1, 1):
        cx, cy = x + d * span, y - 1.05
        ax.plot([x, cx], [y - 0.17, cy + 0.17], color=INK, lw=1.4)
        grow(cx, cy, span / 2, depth + 1, (cx - 0.4) / 6.3)


grow(4.0, 4.4, 1.7, 0, 0.5)
save(fig, "fate.png")

# ---------------- 3. transfer: matching regulatory elements from mouse to human ----------------
fig, ax = canvas()
for y, lab, col in ((3.6, "mouse", "#D98E3A"), (1.4, "human", TEAL)):
    ax.plot([1.35, 7.4], [y, y], color=INK, lw=1.6)
    ax.text(1.25, y, lab, ha="right", va="center", fontsize=20, color=col, fontweight="bold")
    ax.add_patch(FancyBboxPatch((6.3, y - 0.18), 0.9, 0.36, boxstyle="round,pad=0.02", fc=LIGHT, ec=INK, lw=1.2))
    ax.text(6.75, y, "gene", ha="center", va="center", fontsize=16, color=INK)
mouse_el = [2.0, 3.3, 4.7]
human_el = [2.3, 3.8, 5.1]
for x in mouse_el:
    ax.add_patch(Rectangle((x - 0.22, 3.42), 0.44, 0.36, color="#D98E3A", lw=0))
for i, x in enumerate(human_el):
    ax.add_patch(Rectangle((x - 0.22, 1.22), 0.44, 0.36, color=TEAL if i == 1 else LIGHT, lw=0))
# the confident match is drawn solid, the others faint, to say that transfer is a judgment
for i, (mx, hx) in enumerate(zip(mouse_el, human_el)):
    ax.plot([mx, hx], [3.4, 1.6], color=TEAL if i == 1 else GREY, lw=2.2 if i == 1 else 1.0,
            ls="-" if i == 1 else (0, (3, 3)))
ax.text(3.7, 2.5, "same role?", ha="left", va="center", fontsize=18, color=GREY, style="italic")
save(fig, "transfer.png")
print("wrote", sorted(os.listdir(OUT)))
