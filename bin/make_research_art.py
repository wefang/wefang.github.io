"""Draw the three research-area images for the about page, wide and short so each spans the page.

  lineage.png   a lineage tree with a cell at each node; the recorder integrations inside each nucleus
                accumulate edits division by division (colored as in the research statement)
  fate.png      a small simulation: clones grow and divide across a tissue with a signal gradient,
                shown once colored by clone and once colored by the fate each cell took
  transfer.png  a phylogeny with a silhouette per species, a functional genomics signal track for each
                with orthologous genes marked, each organism's proposal onto a human element, and the
                cell states of each species as clusters in a stack of embedding planes

Schematic throughout, and none of it reproduces a figure from a proposal or an unpublished manuscript.

    python bin/make_research_art.py
"""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "research")
os.makedirs(OUT, exist_ok=True)

# the research statement's palette
INK, MUTED, GREY, HAIR = "#1A171B", "#5F5A63", "#BDB7BD", "#E2DDE2"
TEAL, ORANGE, PURPLE, BROWN, MAGENTA, BLUE = "#1B9E77", "#D95F02", "#7570B3", "#A6761D", "#E7298A", "#1F9BD6"
ALLELES = [TEAL, ORANGE, PURPLE, BROWN, MAGENTA, BLUE]
DROPOUT = "#C9C5CA"
plt.rcParams["font.family"] = "DejaVu Sans"
W, H = 12.0, 3.6
DPI = 200


def canvas():
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=DPI, facecolor="white")
    plt.close(fig)


def text(ax, x, y, s, **kw):
    kw.setdefault("fontsize", 12)
    kw.setdefault("color", MUTED)
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, s, **kw)


# =============================== 1. lineage ===============================
# A lineage tree drawn with a cell at every node. Each cell carries the same set of recorder
# integrations inside its nucleus; every division writes a new edit into one integration of a daughter,
# so the edits accumulate down the tree and each sampled cell ends with its own record.
fig, ax = canvas()
rng = np.random.default_rng(9)
NSITE = 5
LEVELS = [3.1, 2.25, 1.4, 0.52]                    # cell centers, root at the top
LEAF_X = np.linspace(1.35, 9.85, 8)
UNEDITED = "#2B2530"


def cell(x, y, sites, rx=0.42, ry=0.33, seed=0):
    """A cell with an irregular outline, and a nucleus holding the integration sites."""
    r = np.random.default_rng(seed)
    th = np.linspace(0, 2 * np.pi, 120)
    wob = 1 + 0.07 * np.sin(3 * th + r.uniform(0, 6)) + 0.05 * np.cos(5 * th + r.uniform(0, 6)) \
        + 0.03 * np.sin(7 * th + r.uniform(0, 6))
    ax.add_patch(Polygon(np.c_[x + rx * wob * np.cos(th), y + ry * wob * np.sin(th)], closed=True,
                         fc="#F6EFF3", ec="#8E8395", lw=1.2, zorder=3))
    nx, ny = x + r.normal(0, 0.02), y + r.normal(0, 0.015)
    ax.add_patch(matplotlib.patches.Ellipse((nx, ny), 0.56, 0.36, angle=r.uniform(-8, 8),
                                            fc="#E6DDEA", ec="#A497AD", lw=0.9, zorder=4))
    w, gap = 0.064, 0.024
    x0 = nx - (NSITE * w + (NSITE - 1) * gap) / 2
    for s in range(NSITE):
        c = sites[s] if sites[s] is not None else UNEDITED
        ax.add_patch(Rectangle((x0 + s * (w + gap), ny - 0.1), w, 0.2, fc=c, ec="none", zorder=5))


def grow(level, lo, hi, sites, px=None, py=None, node=0):
    x = float(np.mean(LEAF_X[lo:hi]))
    y = LEVELS[level]
    if px is not None:                           # the branch from the parent, with the edit it wrote
        ax.plot([px, px, x, x], [py - 0.34, (py + y) / 2 + 0.05, (py + y) / 2 + 0.05, y + 0.34],
                color=INK, lw=1.4, solid_joinstyle="round", zorder=2)
    cell(x, y, sites, seed=node)
    if level == len(LEVELS) - 1:
        return
    mid = (lo + hi) // 2
    for k, (a, b) in enumerate(((lo, mid), (mid, hi))):
        child = list(sites)
        free = [s for s in range(NSITE) if child[s] is None]
        if free and rng.random() < 0.9:
            s = int(rng.choice(free))
            child[s] = ALLELES[int(rng.integers(len(ALLELES)))]
        grow(level + 1, a, b, child, x, y, node * 2 + 1 + k)


grow(0, 0, 8, [None] * NSITE)
ax.add_patch(FancyArrowPatch((0.45, 3.35), (0.45, 0.3), arrowstyle="-|>", mutation_scale=12, color=GREY, lw=1.2))
text(ax, 0.22, 1.85, "Time", fontsize=11, rotation=90)
# key
KX = 10.55
ax.add_patch(Rectangle((KX, 2.62), 0.07, 0.2, fc=UNEDITED, ec="none"))
text(ax, KX + 0.17, 2.72, "Unedited", fontsize=10.5, ha="left")
for j, c in enumerate(ALLELES[:4]):
    ax.add_patch(Rectangle((KX + j * 0.09, 2.2), 0.07, 0.2, fc=c, ec="none"))
text(ax, KX + 0.43, 2.3, "Edited", fontsize=10.5, ha="left")
text(ax, KX + 0.55, 1.8, "Integrations\ninside each\nnucleus", fontsize=9.5, color=MUTED)
save(fig, "lineage.png")

# =============================== 2. fate ===============================
# A growing tissue in physical units (cell diameter 1). Every generation each cell divides, its
# daughter is placed beside it along a slightly preferred axis, and overlapping cells push apart, so
# the tissue expands and each clone stays a coherent, ragged patch. Each cell is drawn as its Voronoi
# tile cut to a disk around the cell, which gives the scalloped edge of a real epithelium, with a
# nucleus. The signal is a gradient along the long axis; at the end each cell reads its fate from the
# signal it sits in, plus a bias inherited from the founder of its clone.
from scipy.spatial import Voronoi
from matplotlib.patches import Polygon as Poly, Ellipse
import matplotlib.image as mpimg
from shapely.geometry import Point, Polygon as SPoly
from shapely.ops import unary_union

rs = np.random.default_rng(21)
GENS, CLONE_GEN, EARLY_GEN = 9, 3, 6


def relax(p, iters=30):
    for _ in range(iters):
        d = p[:, None, :] - p[None, :, :]
        r = np.hypot(d[..., 0], d[..., 1]) + np.eye(len(p)) * 9
        over = np.clip(1.0 - r, 0, None)
        p = p + 0.25 * (d / r[..., None] * over[..., None]).sum(1)
    return p


pos = np.zeros((1, 2))
clone = np.array([-1])
bias = np.zeros(1)
snap = {}
for g in range(1, GENS + 1):
    ang = rs.normal(0, 0.38, len(pos)) + np.where(rs.random(len(pos)) < 0.5, 0, np.pi)   # mostly along x
    off = 0.5 * np.c_[np.cos(ang), 0.7 * np.sin(ang)]
    pos = np.vstack([pos - off, pos + off])
    clone = np.concatenate([clone, clone])
    bias = np.concatenate([bias, bias])
    if g == CLONE_GEN:                                   # label the founders of the clones we follow
        clone = np.arange(len(pos))
        bias = rs.normal(0, 0.07, len(pos))
    pos = relax(pos + rs.normal(0, 0.08, pos.shape))
    if g == EARLY_GEN:
        snap = dict(pos=pos.copy(), clone=clone.copy())

XMIN, XMAX = pos[:, 0].min(), pos[:, 0].max()


def signal_of(p):
    return np.clip((p[:, 0] - XMIN) / (XMAX - XMIN), 0, 1)


fate = np.digitize(signal_of(pos) + bias + rs.normal(0, 0.05, len(pos)), [0.38, 0.64])

CLONE_COLS = ["#8DD3C7", "#FDB462", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#D9C27A"]
FATE_COLS = ["#F2C14E", "#3FA7A0", "#3D2C6B"]
SIG_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("sig", ["#F7F4FB", "#9E9AC8", "#3F007D"])


def tissue(ax, cx, cy, scale, p, colors):
    """Cells as Voronoi tiles cut to a disk around each cell; returns the tissue polygon."""
    xy = np.c_[cx + scale * p[:, 0], cy + scale * p[:, 1]]
    lo, hi = xy.min(0) - 10, xy.max(0) + 10
    far = np.array([[lo[0], lo[1]], [hi[0], lo[1]], [lo[0], hi[1]], [hi[0], hi[1]]])
    vor = Voronoi(np.vstack([xy, far]))
    shapes = []
    for i in range(len(xy)):
        reg = vor.regions[vor.point_region[i]]
        if -1 in reg or not reg:
            continue
        cell = SPoly(vor.vertices[reg]).intersection(Point(xy[i]).buffer(0.62 * scale, 24))
        if cell.is_empty:
            continue
        shapes.append(cell)
        ax.add_patch(Poly(np.array(cell.exterior.coords), closed=True, fc=colors[i], ec="white", lw=0.7, zorder=3))
        c = cell.centroid
        ax.add_patch(Ellipse((c.x, c.y), 0.34 * scale, 0.28 * scale, angle=rs.uniform(0, 180),
                             fc="#2B2530", ec="none", alpha=0.55, zorder=4))
    body = unary_union([c.buffer(0.08 * scale) for c in shapes]).buffer(-0.06 * scale)   # close the gaps between cells
    for g in getattr(body, "geoms", [body]):
        ax.add_patch(Poly(np.array(g.exterior.coords), closed=True, fc="none", ec=MUTED, lw=1.0, zorder=5))
    return body


fig, ax = canvas()
# early: two small tissues, clones and the signal they sit in
early, eclone = snap["pos"], snap["clone"]
# one scale for both stages, set so the grown tissue fills its panel; the early tissue is true to size
SC = min(3.3 / np.ptp(pos[:, 0]), 2.45 / np.ptp(pos[:, 1]))
tissue(ax, 1.15, 1.9, SC, early, [CLONE_COLS[c % 8] for c in eclone])
tissue(ax, 3.2, 1.9, SC, early, [SIG_CMAP(v) for v in signal_of(early)])
# late: the grown tissue, clones and the fate each cell took
tissue(ax, 6.45, 1.9, SC, pos, [CLONE_COLS[c % 8] for c in clone])
tissue(ax, 10.05, 1.9, SC, pos, [FATE_COLS[f] for f in fate])
ax.add_patch(FancyArrowPatch((4.2, 1.9), (4.7, 1.9), arrowstyle="-|>", mutation_scale=14, color=MUTED, lw=1.4))
# group labels on top, panel labels underneath, so the two never read as one
for x0, x1, lab in ((0.3, 4.05, "Early"), (4.85, 11.85, "Late")):
    ax.plot([x0, x1], [3.32, 3.32], color=GREY, lw=1.0)
    text(ax, (x0 + x1) / 2, 3.47, lab, fontsize=12.5, color=INK, fontweight="bold",
         bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none"))
for x, lab in ((1.15, "Clones"), (3.2, "Signal"), (6.45, "Clones"), (10.05, "Fates")):
    text(ax, x, 0.28, lab, fontsize=11.5, color=MUTED)
# fate key
for k, (c, lab) in enumerate(zip(FATE_COLS, ("Fate 1", "Fate 2", "Fate 3"))):
    ax.add_patch(Rectangle((8.95 + k * 0.85, 0.02), 0.16, 0.12, fc=c, ec="none"))
    text(ax, 9.16 + k * 0.85, 0.08, lab, fontsize=9.5, ha="left")
save(fig, "fate.png")

# =============================== 3. transfer ===============================
# Left: a phylogeny with a silhouette per species; for each, a functional genomics signal track over
# the same region, orthologous genes marked, and the region's synteny drawn as conserved blocks along
# the genome (their order and orientation differ between species); each model organism's proposal
# points onto the human track.
# Right: the cell states of the same species as a stack of embedding planes, as in the research
# statement: a branching trajectory in each, corresponding states linked, and human sampled only to
# an earlier stage, so its later states are outlines.
SIL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "silhouettes")
W3, H3 = 12.0, 4.6
fig = plt.figure(figsize=(W3, H3), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W3)
ax.set_ylim(0, H3)
ax.set_aspect("equal")
ax.axis("off")
rs = np.random.default_rng(12)
SPECIES = [("Zebrafish", "zebrafish", BLUE), ("Mouse", "mouse", ORANGE), ("Macaque", "macaque", PURPLE), ("Human", "human", TEAL)]
ROW = [3.45, 2.52, 1.59, 0.66]                   # track baselines, human at the bottom


def hline(x0, x1, y):
    ax.plot([x0, x1], [y, y], color=INK, lw=1.5, solid_capstyle="round")


def vline(x, y0, y1):
    ax.plot([x, x], [y0, y1], color=INK, lw=1.5, solid_capstyle="round")


mid = [y + 0.3 for y in ROW]
n3 = (mid[2] + mid[3]) / 2
n2 = (mid[1] + n3) / 2
n1 = (mid[0] + n2) / 2
hline(0.12, 0.35, n1)
vline(0.35, mid[0], n2)
hline(0.35, 0.95, mid[0])
hline(0.35, 0.55, n2)
vline(0.55, mid[1], n3)
hline(0.55, 0.95, mid[1])
hline(0.55, 0.75, n3)
vline(0.75, mid[2], mid[3])
hline(0.75, 0.95, mid[2])
hline(0.75, 0.95, mid[3])


def silhouette(name, x, y, hmax, wmax, color):
    img = mpimg.imread(os.path.join(SIL, name + ".png"))
    h, w = img.shape[:2]
    sc = min(hmax / h, wmax / w)
    rgba = np.zeros((h, w, 4))
    rgba[..., :3] = matplotlib.colors.to_rgb(color)
    rgba[..., 3] = img[..., 3] if img.shape[2] == 4 else 1 - img[..., 0]
    ax.imshow(rgba, extent=(x - w * sc / 2, x + w * sc / 2, y - h * sc / 2, y + h * sc / 2), zorder=4)


for (lab, key, col), y in zip(SPECIES, ROW):
    silhouette(key, 1.35, y + 0.4, 0.55, 0.7, col)
    text(ax, 1.35, y + 0.02, lab, fontsize=10.5, color=col, fontweight="bold")

TX0, TX1 = 1.95, 6.05
xs = np.linspace(TX0, TX1, 700)
# element positions per species; the region is longer in mouse and shorter in zebrafish, and the
# middle block is inverted in zebrafish, so the elements do not line up column by column
ELEM = {"zebrafish": [2.95, 4.45, 5.0], "mouse": [2.85, 3.95, 5.1], "macaque": [3.0, 4.05, 5.05], "human": [3.02, 4.08, 5.07]}
GA = {"zebrafish": 2.05, "mouse": 2.02, "macaque": 2.08, "human": 2.08}
GB = {"zebrafish": 5.45, "mouse": 5.6, "macaque": 5.55, "human": 5.57}
def gene(x, y, color, label=None, strand=1):
    """A gene arrow; strand -1 points left, as a gene on the reverse strand."""
    pts = [(0, -0.07), (0.26, -0.07), (0.34, 0), (0.26, 0.07), (0, 0.07)]
    if strand < 0:
        pts = [(0.34 - px, py) for px, py in pts]
    ax.add_patch(Polygon([(x + px, y + py) for px, py in pts], closed=True, fc=color, ec="none", zorder=5))
    if label:
        text(ax, x + 0.17, y + 0.2, label, fontsize=9.5, color=INK, style="italic")


for (lab, key, col), y in zip(SPECIES, ROW):
    sig = 0.03 * rs.random(xs.size)
    for k, e in enumerate(ELEM[key]):
        h = 0.5
        sig += h * np.exp(-((xs - e) ** 2) / (2 * 0.045 ** 2))
    for e in rs.uniform(TX0 + 0.5, TX1 - 0.8, 3):                       # species-specific peaks
        sig += 0.15 * np.exp(-((xs - e) ** 2) / (2 * 0.03 ** 2))
    ax.fill_between(xs, y, y + sig, color=col, alpha=0.9, lw=0, zorder=3)
    ax.plot([TX0, TX1], [y, y], color=GREY, lw=1, zorder=2)
    top = key == "zebrafish"
    gene(GA[key], y, INK, "Gene A" if top else None)
    gene(GB[key], y, INK, "Gene B" if top else None, strand=-1)
# orthologous genes joined across species
for (_, a, _), (_, b, _), ya, yb in zip(SPECIES, SPECIES[1:], ROW, ROW[1:]):
    for g in (GA, GB):
        ax.plot([g[a] + 0.17, g[b] + 0.17], [ya - 0.1, yb + 0.1], color=MUTED, lw=0.9, ls=(0, (2, 2)), zorder=1)
# each model organism proposes a different human element, from the element it has evidence for
for (lab, key, col), y, k in zip(SPECIES[:3], ROW[:3], (2, 1, 0)):
    ax.add_patch(FancyArrowPatch((ELEM[key][k] + 0.1, y + 0.3), (ELEM["human"][k] + 0.08, ROW[3] + 0.5),
                                 connectionstyle="arc3,rad=-0.3", arrowstyle="-|>", mutation_scale=11,
                                 color=col, lw=1.4, zorder=6))
text(ax, (TX0 + TX1) / 2, 4.38, "Regulatory elements", fontsize=12, color=INK)

# ---------------- right: a stack of embedding planes, one per species ----------------
# each cell state is a cluster, placed alike in every species with a small species-specific shift;
# matching clusters are linked plane to plane, and the two latest states are missing in human
STATES = [("State 1", "#F28E2B", 0.13, 0.5), ("State 2", "#E15759", 0.38, 0.28), ("State 3", "#76B7B2", 0.38, 0.72),
          ("State 4", "#59A14F", 0.72, 0.25), ("State 5", "#B07AA1", 0.78, 0.72)]


def Pk(k, a, b):
    ox, oy = 6.55 + 0.34 * k, ROW[k] - 0.08
    return ox + 3.5 * a + 0.7 * b, oy + 0.62 * b


CENT = []
for k, (lab, key, col) in enumerate(SPECIES):
    ax.add_patch(Polygon([Pk(k, 0, 0), Pk(k, 1, 0), Pk(k, 1, 1), Pk(k, 0, 1)], closed=True,
                         fc="#FAFAFB", ec=col, lw=1.1, alpha=0.95, zorder=9 - k))
    cents = {}
    for j, (sname, sc, ca, cb) in enumerate(STATES):
        ca2, cb2 = ca + rs.normal(0, 0.025), cb + rs.normal(0, 0.04)
        cents[j] = Pk(k, ca2, cb2)
        if key == "human" and j >= 3:                                   # not sampled in human
            ax.add_patch(Ellipse(cents[j], 0.5, 0.2, fc="none", ec=sc, lw=1.2, ls=(0, (1.5, 1.5)), zorder=11))
            continue
        n = 28
        aa = np.clip(ca2 + rs.normal(0, 0.035, n), 0.03, 0.97)
        bb = np.clip(cb2 + rs.normal(0, 0.09, n), 0.06, 0.94)
        pts = np.array([Pk(k, u, v) for u, v in zip(aa, bb)])
        ax.scatter(pts[:, 0], pts[:, 1], s=8, color=sc, edgecolors="none", zorder=10 - k)
    CENT.append(cents)
for ca, cb in zip(CENT, CENT[1:]):
    for j in range(len(STATES)):
        ax.plot([ca[j][0], cb[j][0]], [ca[j][1], cb[j][1]], color="#9FB6D6", lw=0.9, ls=(0, (3, 2)), zorder=1)
text(ax, 8.9, 4.38, "Cell states, aligned across species", fontsize=12, color=INK)
for j, (sname, sc, _, _) in enumerate(STATES):
    ax.scatter([6.75 + j * 1.0], [0.2], s=26, color=sc, edgecolors="none")
    text(ax, 6.85 + j * 1.0, 0.2, sname, fontsize=9.5, ha="left")
fig.savefig(os.path.join(OUT, "transfer.png"), dpi=DPI, facecolor="white")
plt.close(fig)
print("wrote", sorted(os.listdir(OUT)))
