"""Draw the three research-program images for the about page.

Deliberately schematic: each names the idea of one direction without reproducing any
figure from a proposal or an unpublished manuscript.
  lineage.png   a record of edits accumulates as cells divide, and is read back into a dated tree
  fate.png      sister cells in the same state take different fates at different positions
  transfer.png  cell states and regulatory elements are matched from mouse to human

Palette: greys and one ink-blue accent, to sit on a minimal black-and-white site.

    python bin/make_research_art.py
"""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.path import Path
from matplotlib.patches import PathPatch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "research")
os.makedirs(OUT, exist_ok=True)

INK, G1, G2, G3, G4 = "#1D1D1F", "#6E6E73", "#AEAEB2", "#D8D8DC", "#F2F2F4"
ACC, ACC_L = "#2F4F7F", "#9DB0CC"
plt.rcParams["font.family"] = "DejaVu Sans"
W, H = 9.0, 5.0


def canvas():
    fig = plt.figure(figsize=(W, H), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=200, facecolor="white")
    plt.close(fig)


def label(ax, x, y, s, **kw):
    kw.setdefault("fontsize", 13)
    kw.setdefault("color", G1)
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, s, **kw)


def arrow(ax, a, b, color=G1, lw=1.6):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=16, color=color, lw=lw,
                                 shrinkA=0, shrinkB=0))


# =============================== 1. lineage ===============================
fig, ax = canvas()
# a lineage of divisions on the left; each branch may write one edit into a 6-site record
SITES = 6
rng = np.random.default_rng(11)
EDIT_COLS = [ACC, INK, ACC_L, G1]
tree = {}          # node id -> (x, y, parent, edits)
leaves_x = np.linspace(0.55, 4.05, 8)


def build(node, depth, xs, y, parent_edits):
    edits = dict(parent_edits)
    if depth > 0:                                     # this branch writes one new edit
        free = [s for s in range(SITES) if s not in edits]
        if free and rng.random() < 0.85:
            edits[int(rng.choice(free))] = EDIT_COLS[int(rng.integers(len(EDIT_COLS)))]
    x = float(np.mean(xs))
    tree[node] = (x, y, edits)
    if depth < 3:
        h = len(xs) // 2
        build(node * 2 + 1, depth + 1, xs[:h], y - 1.0, edits)
        build(node * 2 + 2, depth + 1, xs[h:], y - 1.0, edits)


build(0, 0, list(leaves_x), 4.35, {})
for node, (x, y, edits) in tree.items():
    for c in (node * 2 + 1, node * 2 + 2):
        if c in tree:
            cx, cy, cedits = tree[c]
            ax.plot([x, cx], [y - 0.13, cy + 0.13], color=G2, lw=1.4, solid_capstyle="round")
            new = set(cedits) - set(edits)
            if new:                                   # the edit, marked where it was written
                mx, my = (x + cx) / 2, (y + cy) / 2
                ax.add_patch(Rectangle((mx - 0.06, my - 0.06), 0.12, 0.12, color=cedits[list(new)[0]], lw=0, zorder=4))
    ax.add_patch(Circle((x, y), 0.13, facecolor="white", edgecolor=INK, lw=1.3, zorder=3))
# the record each sampled cell carries
for node, (x, y, edits) in tree.items():
    if node >= 7:
        for s in range(SITES):
            col = edits.get(s)
            ax.add_patch(Rectangle((x - 0.21 + s * 0.07, 0.62), 0.062, 0.42,
                                   facecolor=col if col else G4, edgecolor="none", zorder=3))
label(ax, 2.3, 0.3, "DNA record in each cell", fontsize=12)

arrow(ax, (4.3, 2.6), (4.95, 2.6), color=INK)
label(ax, 4.62, 2.95, "infer", fontsize=12)

# the dated tree read back from the records
T0, T1 = 4.35, 0.85
xs = np.linspace(5.9, 8.6, 8)
heights = {0: 4.35, 1: 3.55, 2: 3.15, 3: 2.45, 4: 1.95, 5: 2.75, 6: 2.2}
def dated(node, lo, hi):
    if node >= 7:
        x = float(np.mean(xs[lo:hi]))
        ax.add_patch(Circle((x, T1), 0.09, facecolor="white", edgecolor=INK, lw=1.2, zorder=3))
        return x, T1
    mid = (lo + hi) // 2
    xl, yl = dated(node * 2 + 1, lo, mid)
    xr, yr = dated(node * 2 + 2, mid, hi)
    y = heights[node]
    ax.plot([xl, xl, xr, xr], [yl, y, y, yr], color=INK, lw=1.5, solid_joinstyle="round")
    ax.add_patch(Circle(((xl + xr) / 2, y), 0.05, color=ACC, zorder=3))
    return (xl + xr) / 2, y
dated(0, 0, 8)
# time axis, left of the tree, reading down from the first division
ax.add_patch(FancyArrowPatch((5.45, T0), (5.45, T1), arrowstyle="-|>", mutation_scale=12, color=G2, lw=1.1))
label(ax, 5.25, (T0 + T1) / 2, "time", fontsize=11, rotation=90)
label(ax, 7.25, 0.3, "dated lineage tree", fontsize=12)
save(fig, "lineage.png")

# =============================== 2. fate ===============================
fig, ax = canvas()
# the signal: a field across the bottom, pale on the left and dense on the right
field = np.linspace(0, 1, 400).reshape(1, -1)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("sig", ["#FFFFFF", ACC_L, ACC])
ax.imshow(field, extent=(0.5, 8.5, 0.55, 0.9), cmap=cmap, aspect="auto", zorder=1)
label(ax, 0.5, 0.3, "low signal", fontsize=11, ha="left")
label(ax, 8.5, 0.3, "high signal", fontsize=11, ha="right")

FATE = [G4, G2, INK]                  # three fates, by where a cell ends up
def fate_of(x):
    f = (x - 0.5) / 8.0
    return 0 if f < 0.36 else (1 if f < 0.66 else 2)

cells = {}
def grow(node, x, y, span, depth):
    cells[node] = (x, y, depth)
    if depth < 3:
        for i, d in enumerate((-1, 1)):
            grow(node * 2 + 1 + i, x + d * span, y - 1.05, span * 0.5, depth + 1)
grow(0, 4.5, 4.4, 2.0, 0)
for node, (x, y, depth) in cells.items():
    for c in (node * 2 + 1, node * 2 + 2):
        if c in cells:
            cx, cy, _ = cells[c]
            verts = [(x, y - 0.17), (x, (y + cy) / 2), (cx, (y + cy) / 2 + 0.25), (cx, cy + 0.17)]
            ax.add_patch(PathPatch(Path(verts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                   fc="none", ec=G2, lw=1.4))
for node, (x, y, depth) in cells.items():
    if depth == 3:
        ax.add_patch(Circle((x, y), 0.19, facecolor=FATE[fate_of(x)], edgecolor=INK, lw=1.2, zorder=3))
        ax.plot([x, x], [y - 0.19, 0.92], color=G3, lw=0.9, ls=(0, (2, 2)), zorder=0)
    else:
        ax.add_patch(Circle((x, y), 0.19, facecolor="white", edgecolor=INK, lw=1.2, zorder=3))
# two sisters, one state, two fates: the question the model is built to answer
a, b = cells[9], cells[10]
ax.add_patch(FancyBboxPatch((a[0] - 0.42, a[1] - 0.38), (b[0] - a[0]) + 0.84, 0.76,
                            boxstyle="round,pad=0.02,rounding_size=0.18", fc="none", ec=ACC, lw=1.4, ls=(0, (4, 3))))
label(ax, (a[0] + b[0]) / 2, a[1] + 0.62, "same state, different fate", fontsize=11.5, color=ACC,
      bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"), zorder=5)
save(fig, "fate.png")

# =============================== 3. transfer ===============================
fig, ax = canvas()


def trajectory(ox, oy, human):
    """A small branching path of cell states. The human one is sampled to an earlier stage,
    so its branches stop short and continue as an outline."""
    stem = [(ox, oy), (ox + 0.9, oy)]
    ax.plot(*zip(*stem), color=INK, lw=1.8, solid_capstyle="round")
    for dy, stop in ((0.62, 0.55), (-0.62, 0.55)):
        t = np.linspace(0, 1, 40)
        xs = ox + 0.9 + 1.9 * t
        ys = oy + dy * (3 * t ** 2 - 2 * t ** 3)
        cut = 40 if not human else int(40 * stop)
        ax.plot(xs[:cut], ys[:cut], color=INK, lw=1.8, solid_capstyle="round")
        if human:
            ax.plot(xs[cut - 1:], ys[cut - 1:], color=G2, lw=1.5, ls=(0, (2, 2.5)))
        ax.add_patch(Circle((xs[-1], ys[-1]), 0.12, facecolor="white" if human else INK, edgecolor=INK if not human else G1, lw=1.2, zorder=3))
    ax.add_patch(Circle((ox + 0.9, oy), 0.12, facecolor=INK, edgecolor=INK, zorder=3))


def genome(ox, oy, human):
    ax.plot([ox, ox + 3.2], [oy, oy], color=INK, lw=1.5)
    els = [ox + 0.45, ox + 1.2, ox + 1.95]
    for i, x in enumerate(els):
        fill = (ACC if i == 1 else G3) if human else INK
        ax.add_patch(Rectangle((x - 0.17, oy - 0.15), 0.34, 0.30, facecolor=fill, edgecolor="none", zorder=3))
    ax.add_patch(FancyBboxPatch((ox + 2.55, oy - 0.16), 0.6, 0.32, boxstyle="round,pad=0.02,rounding_size=0.06",
                                fc="white", ec=INK, lw=1.2, zorder=3))
    label(ax, ox + 2.85, oy, "gene", fontsize=9.5, color=INK)
    return els


label(ax, 2.75, 4.55, "mouse", fontsize=14, color=INK, fontweight="bold")
label(ax, 7.05, 4.55, "human", fontsize=14, color=INK, fontweight="bold")
trajectory(1.15, 3.3, human=False)
trajectory(5.45, 3.3, human=True)
m_el = genome(1.15, 1.35, human=False)
h_el = genome(5.45, 1.35, human=True)
label(ax, 0.5, 3.3, "cell\nstates", fontsize=10, color=G1)
label(ax, 0.5, 1.35, "regulatory\nelements", fontsize=10, color=G1)
# matched states (dashed) and matched elements, the confident match drawn solid in the accent
for dy in (0.62, -0.62):
    ax.add_patch(FancyArrowPatch((1.15 + 2.8 + 0.14, 3.3 + dy), (5.45 + 2.8 - 0.14, 3.3 + dy),
                                 connectionstyle="arc3,rad=-0.10", arrowstyle="-", color=G2, lw=1.1, ls=(0, (3, 3))))
for i, (mx, hx) in enumerate(zip(m_el, h_el)):
    ax.add_patch(FancyArrowPatch((mx, 1.52), (hx, 1.52), connectionstyle="arc3,rad=-0.28", arrowstyle="-",
                                 color=ACC if i == 1 else G3, lw=2.0 if i == 1 else 1.0,
                                 ls="-" if i == 1 else (0, (3, 3))))
label(ax, 4.9, 0.55, "which human element plays the same role?", fontsize=11.5, color=ACC)
save(fig, "transfer.png")
print("wrote", sorted(os.listdir(OUT)))
