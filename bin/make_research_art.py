"""Draw the three research-area images for the about page, wide and short so each spans the page.

  lineage.png   a lineage tree with a cell at each node; the recorder integrations inside each nucleus
                accumulate edits division by division, the newest one ringed, and the sampled cells are
                sequenced into a character matrix (colored as in the research statement)
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

# ---- one palette for all three figures ----
# Categorical colors are the Okabe-Ito set, which stays distinguishable under common color vision
# deficiencies. The same meaning keeps the same color across figures: cell identities (fates in the
# fate figure, states in the transfer figure) share IDENT, and recorder edits reuse its hues. Lightness
# separates layers of meaning: clones are pale tints, identities are saturated. Magnitude (the signal)
# is a single-hue slate ramp, never a category color. One accent, vermillion, marks what matters:
# the newest edit, and human. Everything else is neutral.
INK, MUTED, GREY, HAIR = "#1F2328", "#6B7280", "#B8BDC4", "#E5E7EB"
OI_ORANGE, OI_SKY, OI_GREEN, OI_BLUE, OI_PURPLE = "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7"
OI_YELLOW, OI_VERMILLION = "#F0E442", "#D55E00"
IDENT = [OI_ORANGE, OI_SKY, OI_GREEN, OI_BLUE, OI_PURPLE]
ACCENT = OI_VERMILLION
SLATE = "#3B4A63"
ALLELES = IDENT
DROPOUT = "#CDD1D6"


def mute(c, f=0.38, grey="#A7ABB2"):
    """c pulled toward a warm grey by f, for a quieter version of the same hue."""
    a, g = np.array(matplotlib.colors.to_rgb(c)), np.array(matplotlib.colors.to_rgb(grey))
    return tuple(a + (g - a) * f)


def shade(c, f=0.3):
    """c darkened by f, for outlines and nuclei drawn in the cell's own hue."""
    return tuple(np.array(matplotlib.colors.to_rgb(c)) * (1 - f))


def tint(c, f=0.55):
    """c mixed toward white by f, for a lighter layer of the same hue."""
    a = np.array(matplotlib.colors.to_rgb(c))
    return tuple(a + (1 - a) * f)


# cell identities (fates, and cell states in the transfer figure): the blue, green and red of the
# research statement's Figure 2, then two more muted hues for states 4 and 5
IDENT_M = ["#6F8CA8", "#7DA58F", "#C48B66", "#D2AE5E", "#8F82AE"]
# clones: four hue families that sit away from the fates' blue, green and red and the signal's slate
# (gold, violet, rose, grey), each in a light and a dark tone. Sister clones get neighbouring indices
# and grow side by side, so the order alternates hue and lightness to keep every border visible.
CLONE_SET = ["#E3B94F", "#6E5A9E", "#E39AB8", "#C9C8CF", "#F3DF9B", "#9C4F78", "#B7A6DA", "#6F6D78"]
plt.rcParams["font.family"] = "DejaVu Sans"
W, H = 12.0, 4.1                                   # every figure: 12 wide, with a header band on top
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
    if name == "fate.png" and os.environ.get("FATE_OUT"):
        name = os.environ["FATE_OUT"]
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    crop_to_ink(path)


def crop_to_ink(path, pad=24):
    """Crop a saved figure to its drawn pixels plus one fixed margin, so every figure's content spans
    the full width it is shown at."""
    from PIL import Image
    im = Image.open(path).convert("RGB")
    a = np.asarray(im).astype(int)
    ink = np.where((a < 245).any(axis=2))
    y0, y1, x0, x1 = ink[0].min(), ink[0].max(), ink[1].min(), ink[1].max()
    im.crop((max(0, x0 - pad), max(0, y0 - pad), min(im.width, x1 + pad), min(im.height, y1 + pad))).save(path)


def badge(ax, x, y, n, r=0.13):
    """A small numbered disc, the same for a sampled cell and its row of the sequenced records."""
    ax.add_patch(matplotlib.patches.Circle((x, y), r, fc="#EEF0F3", ec=GREY, lw=0.8, zorder=6))
    ax.text(x, y, str(n), fontsize=9, color=INK, ha="center", va="center", fontweight="bold", zorder=7)


def text(ax, x, y, s, **kw):
    kw.setdefault("fontsize", 12)
    kw.setdefault("color", MUTED)
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, s, **kw)


def header(ax, x0, x1, y, label):
    """The one header style all three figures share: a bold title over a thin rule."""
    ax.plot([x0, x1], [y - 0.19, y - 0.19], color=GREY, lw=1.0)
    text(ax, (x0 + x1) / 2, y, label, fontsize=12.5, color=INK, fontweight="bold")


# =============================== 1. lineage ===============================
# Left: a lineage tree drawn with a cell at every division and every sampled
# tip, division times varying and the tree unbalanced, as a real one is. Each cell carries the same
# integrations in its nucleus; each division writes a new edit into one integration of a daughter
# (ringed), so the edits accumulate down the tree. Right: the sampled cells are sequenced, and their
# records are read out as the rows of a character matrix, with the dropouts a real readout has.
fig, ax = canvas()
rng = np.random.default_rng(9)
NSITE = 5
UNEDITED = "#2B2F36"
NEW = ACCENT
TIPY = 0.6

# ---------------- the tree ----------------
TREE = {
    "r": (4.33, 3.2, ["A", "B"]),
    "A": (2.68, 2.3, ["A1", "t3"]),
    "A1": (2.01, 1.45, ["t1", "t2"]),
    "B": (6.10, 2.42, ["B1", "B2"]),
    "B1": (5.07, 1.55, ["t4", "t5"]),
    "B2": (7.09, 1.45, ["t6", "t7"]),
}
TIPS = ["t1", "t2", "t3", "t4", "t5", "t6", "t7"]
TIPX = dict(zip(TIPS, [1.40, 2.41, 3.47, 4.57, 5.58, 6.58, 7.59]))


def cell(x, y, sites, new=None, rx=0.36, ry=0.29, seed=0):
    """A cell with an irregular outline, and a nucleus holding the integration sites."""
    r = np.random.default_rng(seed)
    th = np.linspace(0, 2 * np.pi, 120)
    wob = 1 + 0.07 * np.sin(3 * th + r.uniform(0, 6)) + 0.05 * np.cos(5 * th + r.uniform(0, 6)) \
        + 0.03 * np.sin(7 * th + r.uniform(0, 6))
    ax.add_patch(Polygon(np.c_[x + rx * wob * np.cos(th), y + ry * wob * np.sin(th)], closed=True,
                         fc="#F6EFF3", ec="#8E8395", lw=1.1, zorder=3))
    nx, ny = x + r.normal(0, 0.02), y + r.normal(0, 0.015)
    ax.add_patch(matplotlib.patches.Ellipse((nx, ny), 0.5, 0.31, angle=r.uniform(-8, 8),
                                            fc="#E6DDEA", ec="#A497AD", lw=0.8, zorder=4))
    w, gap = 0.056, 0.022
    x0 = nx - (NSITE * w + (NSITE - 1) * gap) / 2
    for s in range(NSITE):
        c = sites[s] if sites[s] is not None else UNEDITED
        ax.add_patch(Rectangle((x0 + s * (w + gap), ny - 0.085), w, 0.17, fc=c, ec="none", zorder=5))
        if s == new:                                            # the edit this cell's division wrote
            ax.add_patch(Rectangle((x0 + s * (w + gap) - 0.022, ny - 0.11), w + 0.044, 0.22, fc="none",
                                   ec=NEW, lw=1.6, zorder=6))


def branch(x0, y0, x1, y1):
    t = np.linspace(0, 1, 40)
    ya, yb = y0 - 0.28, y1 + 0.28
    ax.plot(x0 + (x1 - x0) * (3 * t ** 2 - 2 * t ** 3), ya + (yb - ya) * t, color=INK, lw=1.3,
            solid_capstyle="round", zorder=2)


RECORD = {}


def grow(name, sites, new, seed):
    RECORD[name] = sites
    x, y, kids = TREE[name] if name in TREE else (TIPX[name], TIPY, [])
    cell(x, y, sites, new, seed=seed)
    for k, kid in enumerate(kids):
        child = list(sites)
        free = [s for s in range(NSITE) if child[s] is None]
        knew = None
        if free:
            knew = int(rng.choice(free))
            child[knew] = ALLELES[int(rng.integers(len(ALLELES)))]
        kx, ky = (TREE[kid][0], TREE[kid][1]) if kid in TREE else (TIPX[kid], TIPY)
        branch(x, y, kx, ky)
        grow(kid, child, knew, seed * 2 + 1 + k)


grow("r", [None] * NSITE, None, 1)
for i, t in enumerate(TIPS):
    badge(ax, TIPX[t], TIPY - 0.45, i + 1)
ax.add_patch(FancyArrowPatch((0.45, 3.4), (0.45, 0.35), arrowstyle="-|>", mutation_scale=12, color=GREY, lw=1.2))
text(ax, 0.27, 1.9, "Time", fontsize=10.5, rotation=90)
# the tree's key sits with the tree, in the space beside the root
KX, KY = 0.75, 3.28
ax.add_patch(Rectangle((KX, KY - 0.08), 0.06, 0.16, fc=UNEDITED, ec="none"))
text(ax, KX + 0.14, KY, "Unedited", fontsize=9.5, ha="left")
for j, c in enumerate(ALLELES[:3]):
    ax.add_patch(Rectangle((KX + j * 0.075, KY - 0.43), 0.06, 0.16, fc=c, ec="none"))
text(ax, KX + 0.29, KY - 0.35, "Edited", fontsize=9.5, ha="left")
ax.add_patch(Rectangle((KX, KY - 0.78), 0.06, 0.16, fc=IDENT[0], ec="none"))
ax.add_patch(Rectangle((KX - 0.025, KY - 0.805), 0.11, 0.21, fc="none", ec=NEW, lw=1.6))
text(ax, KX + 0.14, KY - 0.7, "Newest edit", fontsize=9.5, ha="left")

# ---------------- sequencing ----------------
ax.add_patch(FancyArrowPatch((8.85, 1.85), (9.45, 1.85), arrowstyle="-|>", mutation_scale=15, color=INK, lw=1.5))
text(ax, 9.15, 2.13, "Sequence", fontsize=10.5, color=INK)
MX, MY, cw, ch = 10.0, 2.95, 0.3, 0.29
for i, t in enumerate(TIPS):
    y = MY - i * (ch + 0.06)
    badge(ax, MX - 0.22, y + ch / 2, i + 1)
    for s in range(NSITE):
        c = RECORD[t][s] if RECORD[t][s] is not None else UNEDITED
        if rng.random() < 0.2:
            c = DROPOUT                                        # the readout failed at this site
        ax.add_patch(Rectangle((MX + s * (cw + 0.05), y), cw, ch, fc=c, ec="none", zorder=3))
header(ax, 0.3, 8.6, H - 0.18, "Lineage")
header(ax, 9.7, 11.85, H - 0.18, "Sequenced records")
ax.add_patch(Rectangle((MX, 0.2), 0.16, 0.16, fc=DROPOUT, ec="none"))
text(ax, MX + 0.26, 0.28, "Dropout", fontsize=9.5, ha="left")
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


# the tissue grows inside an outline that keeps its shape and scales with the number of cells
TISSUE_SHAPE = os.environ.get("TISSUE_SHAPE", "bud")


def outline_unit(kind):
    """A closed outline of unit area, centred on the origin."""
    th = np.linspace(0, 2 * np.pi, 240, endpoint=False)
    if kind == "bean":                                   # kidney: long, with one concave side
        x = 1.55 * np.cos(th)
        y = 0.8 * np.sin(th) - 0.28 * np.cos(2 * th) * (np.sin(th) > 0) - 0.1 * np.cos(th) ** 2
    elif kind == "bud":                                  # limb bud: a broad base and a rounded, uneven outgrowth
        wob = 1 + 0.07 * np.sin(3 * th + 0.7) + 0.045 * np.cos(5 * th + 0.3) + 0.03 * np.sin(7 * th + 1.9)
        x = 1.15 * np.cos(th) * (1 - 0.2 * np.sin(th)) * wob
        dome = (1.45 * np.sin(th) + 0.22 * np.cos(th) * np.sin(th)) * wob - 0.32     # leans to one side
        base = -0.32 + 0.06 * np.sin(th) + 0.035 * np.sin(4 * th + 0.5)               # not quite flat
        y = np.where(np.sin(th) < 0, base, dome)
    else:                                                # lobed, an irregular organic outline
        r = 1 + 0.16 * np.sin(3 * th + 0.4) + 0.1 * np.cos(5 * th + 1.1) + 0.06 * np.sin(2 * th)
        x, y = 1.4 * r * np.cos(th), 0.85 * r * np.sin(th)
    poly = SPoly(np.c_[x, y]).buffer(0.04).buffer(-0.04)          # smooth out any self-crossing
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    c = poly.centroid
    xy = np.array(poly.exterior.coords) - [c.x, c.y]
    return xy / np.sqrt(SPoly(xy).area)


UNIT = outline_unit(TISSUE_SHAPE)


def confine(p):
    """Pull any cell that has left the outline back just inside it; the outline is scaled so its
    area matches the number of cells."""
    shape = SPoly(UNIT * np.sqrt(0.8 * len(p)))       # slightly under the packing area, so cells fill it
    inner = shape.buffer(-0.3)
    for i, q in enumerate(p):
        pt = Point(q)
        if not inner.contains(pt):
            b = inner.exterior.interpolate(inner.exterior.project(pt))
            p[i] = [b.x, b.y]
    return p


def relax(p, iters=30):
    """Overlapping cells push apart, inside the tissue outline."""
    for it in range(iters):
        d = p[:, None, :] - p[None, :, :]
        r = np.hypot(d[..., 0], d[..., 1]) + np.eye(len(p)) * 9
        over = np.clip(1.0 - r, 0, None)
        p = p + 0.45 * (d / r[..., None] * over[..., None]).sum(1)
        if len(p) >= 8 and it % 2 == 1:
            p = confine(p)
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
    pos = relax(pos + rs.normal(0, 0.08, pos.shape), iters=60 + 25 * g)   # enough to resolve the overlaps as the tissue grows
    if g == EARLY_GEN:
        snap = dict(pos=pos.copy(), clone=clone.copy())

# a gentle bend, so the tissue is not a plain ellipse
XMIN, XMAX = pos[:, 0].min(), pos[:, 0].max()


def signal_of(p):
    """Position along the tissue mapped to signal: it changes fast near both ends and hardly at all
    through the middle, where cells sit on a plateau."""
    u = np.clip((p[:, 0] - p[:, 0].min()) / np.ptp(p[:, 0]), 0, 1)       # the gradient spans the tissue it is in
    v = 2 * u - 1
    return 0.5 + 0.5 * np.sign(v) * np.abs(v) ** 3


fate = np.digitize(signal_of(pos) + 0.5 * bias + rs.normal(0, 0.025, len(pos)), [0.465, 0.535])   # the plateau stays uncommitted

CLONE_COLS = CLONE_SET
# low signal, middle, high signal: the middle band stays uncommitted (blue), flanked by two fates
FATE_COLS = [IDENT_M[1], IDENT_M[0], IDENT_M[2]]
FATE_LABS = ("Fate A", "Uncommitted", "Fate B")
# the signal is a single slate ramp, from near white to near black, so its magnitude reads at a glance
SIG_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("sig", ["#F7F8FA", "#B4BDCB", "#5E6B82", "#18202E"])


def tissue(ax, cx, cy, scale, p, colors):
    """Flat 2D cells without walls: each cell's Voronoi tile, cut to a disk around it, touches its
    neighbours, and a flat nucleus marks each cell; the tissue behind them follows the outer cells."""
    xy = np.c_[cx + scale * p[:, 0], cy + scale * p[:, 1]]
    lo, hi = xy.min(0) - 10, xy.max(0) + 10
    far = np.array([[lo[0], lo[1]], [hi[0], lo[1]], [lo[0], hi[1]], [hi[0], hi[1]]])
    vor = Voronoi(np.vstack([xy, far]))
    cells_, keep = [], []
    for i in range(len(xy)):
        reg = vor.regions[vor.point_region[i]]
        if -1 in reg or not reg:
            continue
        cells_.append(SPoly(vor.vertices[reg]).intersection(Point(xy[i]).buffer(0.62 * scale, 24)))
        keep.append(i)
    body = unary_union([c.buffer(0.3 * scale) for c in cells_]).buffer(-0.18 * scale)   # a smooth tissue edge
    for g in getattr(body, "geoms", [body]):
        ax.add_patch(Poly(np.array(g.exterior.coords), closed=True, fc="#F1EEEA", ec="#B5AFA7", lw=1.0, zorder=2))
    for i, c in zip(keep, cells_):
        for g in getattr(c, "geoms", [c]):                        # no wall: neighbours touch
            ax.add_patch(Poly(np.array(g.exterior.coords), closed=True, fc=colors[i], ec=colors[i], lw=0.4, zorder=3))
        ctr = c.centroid
        ax.add_patch(matplotlib.patches.Circle((ctr.x, ctr.y), 0.17 * scale, fc=shade(colors[i], 0.3), ec="none", zorder=4))
    return body


fig, ax = canvas()
# early: two small tissues, clones and the signal they sit in
early, eclone = snap["pos"], snap["clone"]
# one scale for both stages, set so the grown tissue fills its panel; the early tissue is true to size
SC = min(3.35 / np.ptp(pos[:, 0]), 2.75 / np.ptp(pos[:, 1]))
tissue(ax, 1.15, 1.9, SC, early, [CLONE_COLS[c % 8] for c in eclone])
tissue(ax, 3.2, 1.9, SC, early, [SIG_CMAP(v) for v in signal_of(early)])
# late: the grown tissue, clones and the fate each cell took
tissue(ax, 6.45, 1.9, SC, pos, [CLONE_COLS[c % 8] for c in clone])
tissue(ax, 10.05, 1.9, SC, pos, [FATE_COLS[f] for f in fate])
ax.add_patch(FancyArrowPatch((4.2, 1.9), (4.7, 1.9), arrowstyle="-|>", mutation_scale=14, color=MUTED, lw=1.4))
# group labels on top, panel labels underneath, so the two never read as one
for x0, x1, lab in ((0.3, 4.05, "Early"), (4.85, 11.85, "Late")):
    header(ax, x0, x1, H - 0.18, lab)
for x, lab in ((1.15, "Clones"), (3.2, "Signal"), (6.45, "Clones"), (10.05, "Fates")):
    text(ax, x, 0.28, lab, fontsize=11.5, color=MUTED)
# the signal's scale bar, under the signal panel
sbar = np.linspace(0, 1, 200).reshape(1, -1)
ax.imshow(sbar, extent=(2.55, 3.85, 0.02, 0.12), aspect="auto", cmap=SIG_CMAP, zorder=3)
text(ax, 2.5, 0.07, "Low", fontsize=9, ha="right")
text(ax, 3.9, 0.07, "High", fontsize=9, ha="left")
# fate key
for k, (c, lab) in enumerate(zip(FATE_COLS, FATE_LABS)):
    kx = (8.55, 9.45, 10.8)[k]
    ax.add_patch(Rectangle((kx, 0.02), 0.16, 0.12, fc=c, ec="none"))
    text(ax, kx + 0.21, 0.08, lab, fontsize=9.5, ha="left")
save(fig, "fate.png")

# =============================== 3. transfer ===============================
# Left: a phylogeny with a silhouette per species; for each, a functional genomics signal track over
# the same region on its own plane, stacked in perspective, with the orthologous genes marked; each
# model organism proposes a different human element.
# Right: the cell states of the same species as a stack of embedding planes, one cluster per state.
# Human is the reference and has every state; each model organism lacks some of them (outlined).
# Model organisms are drawn in neutrals and human in the accent, so the figure reads human-centered.
SIL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "silhouettes")
W3, H3, Y3 = 12.0, 4.6, 0.2                      # room under the human plane for one label
fig = plt.figure(figsize=(W3, H3 - Y3), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W3)
ax.set_ylim(Y3, H3)
ax.set_aspect("equal")
ax.axis("off")
rs = np.random.default_rng(12)
SPECIES = [("Zebrafish", "zebrafish"), ("Mouse", "mouse"), ("Macaque", "macaque"), ("Human", "human")]
ROW = [3.45, 2.52, 1.59, 0.66]                   # plane bottoms, human at the bottom


def scol(key):
    return ACCENT if key == "human" else SLATE


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


for (lab, key), y in zip(SPECIES, ROW):
    silhouette(key, 1.35, y + 0.4, 0.55, 0.7, scol(key))
    text(ax, 1.35, y + 0.02, lab, fontsize=10.5, color=scol(key), fontweight="bold")


def plane(k, a, b, x0, wa, wb, hb):
    """A point on plane k: a runs along the plane, b runs into it (drawn up and to the right). The
    planes of a stack are aligned vertically, so a link between two planes is a vertical line."""
    return x0 + wa * a + wb * b, ROW[k] - 0.1 + hb * b


def Lk(k, a, b):
    return plane(k, a, b, 2.0, 3.85, 0.55, 0.42)


def Rk(k, a, b):
    return plane(k, a, b, 6.85, 3.75, 0.7, 0.62)


def draw_plane(P, k, key, z):
    ax.add_patch(Polygon([P(k, 0, 0), P(k, 1, 0), P(k, 1, 1), P(k, 0, 1)], closed=True, fc="#F7F7F8",
                         ec=scol(key), lw=1.3 if key == "human" else 0.9, alpha=0.96, zorder=z))


# ---------------- left: a signal track on each plane ----------------
ELEM = {"zebrafish": [0.24, 0.6, 0.73], "mouse": [0.22, 0.48, 0.76], "macaque": [0.26, 0.5, 0.75], "human": [0.26, 0.51, 0.755]}
GA = {"zebrafish": 0.03, "mouse": 0.02, "macaque": 0.04, "human": 0.04}      # gene A, left flank
GB = {"zebrafish": 0.86, "mouse": 0.89, "macaque": 0.88, "human": 0.885}     # gene B, reverse strand
TRACK_B = 0.3                                    # the track runs along the plane at this depth
NOMAP = 0.4                                      # a zebrafish element that maps to nothing in human
aa = np.linspace(0.0, 1.0, 600)
for k, (lab, key) in enumerate(SPECIES):
    draw_plane(Lk, k, key, 2 + k * 0.01)
    sig = 0.03 * rs.random(aa.size)
    for e in ELEM[key]:
        sig += 0.5 * np.exp(-((aa - e) ** 2) / (2 * 0.011 ** 2))
    for e in rs.uniform(0.12, 0.8, 3):                                   # species-specific peaks
        if key == "human" and abs(e - NOMAP) < 0.07:
            continue
        sig += 0.15 * np.exp(-((aa - e) ** 2) / (2 * 0.008 ** 2))
    if key == "zebrafish":                                              # an element with no human counterpart
        sig += 0.42 * np.exp(-((aa - NOMAP) ** 2) / (2 * 0.011 ** 2))
    base = np.array([Lk(k, a, TRACK_B) for a in aa])
    ax.fill_between(base[:, 0], base[:, 1], base[:, 1] + 0.95 * sig, color=scol(key), alpha=0.85, lw=0, zorder=5)
    ax.plot(base[:, 0], base[:, 1], color=GREY, lw=0.9, zorder=4)
    for g, strand, name in ((GA, 1, "Gene A"), (GB, -1, "Gene B")):
        gx, gy = Lk(k, g[key], TRACK_B)
        pts = [(0, -0.065), (0.24, -0.065), (0.32, 0), (0.24, 0.065), (0, 0.065)]
        if strand < 0:
            pts = [(0.32 - px, py) for px, py in pts]
        ax.add_patch(Polygon([(gx + px, gy + py) for px, py in pts], closed=True, fc=INK, ec="none", zorder=6))
        if k == 0:
            text(ax, gx + 0.16, gy + 0.2, name, fontsize=9.5, color=INK, style="italic")
# orthologous genes joined plane to plane
for k in range(3):
    for g in (GA, GB):
        p, q = Lk(k, g[SPECIES[k][1]], TRACK_B), Lk(k + 1, g[SPECIES[k + 1][1]], TRACK_B)
        ax.plot([p[0] + 0.16, q[0] + 0.16], [p[1] - 0.07, q[1] + 0.07], color=MUTED, lw=0.9, ls=(0, (2, 2)), zorder=5.5)
# each model organism proposes a different human element, from the element it has evidence for
for k, j in zip(range(3), (2, 1, 0)):
    key = SPECIES[k][1]
    p = Lk(k, ELEM[key][j], TRACK_B)
    q = Lk(3, ELEM["human"][j], TRACK_B)
    ax.add_patch(FancyArrowPatch((p[0] + 0.05, p[1] + 0.3), (q[0] + 0.03, q[1] + 0.5),
                                 connectionstyle="arc3,rad=-0.3", arrowstyle="-|>", mutation_scale=11,
                                 color=INK, lw=1.3, zorder=7))
    hx, hy = q
    ax.add_patch(matplotlib.patches.Ellipse((hx, hy + 0.22), 0.2, 0.56, fc="none", ec=ACCENT, lw=1.4, zorder=7))
# the no-map example: the zebrafish element has no counterpart in human, so its line ends in a cross
p = Lk(0, NOMAP, TRACK_B)
q = Lk(3, NOMAP, TRACK_B)
ax.plot([p[0], q[0] + 0.02], [p[1] - 0.02, q[1] + 0.12], color=MUTED, lw=1.2, ls=(0, (3, 2)), zorder=6)
ax.plot([q[0] - 0.07, q[0] + 0.11], [q[1] + 0.05, q[1] + 0.23], color=MUTED, lw=2.0, zorder=7)
ax.plot([q[0] - 0.07, q[0] + 0.11], [q[1] + 0.23, q[1] + 0.05], color=MUTED, lw=2.0, zorder=7)
text(ax, q[0] - 0.1, ROW[3] - 0.26, "No human counterpart", fontsize=9, color=MUTED)
header(ax, 1.95, 6.65, H3 - 0.18, "Regulatory elements")

# ---------------- right: a stack of embedding planes, one cluster per cell state ----------------
STATES = [("State %d" % (j + 1), IDENT_M[j], ca, cb) for j, (ca, cb) in
          enumerate(((0.13, 0.5), (0.4, 0.3), (0.62, 0.72), (0.85, 0.35)))]
# the missing states follow the phylogeny: the more distant the species, the more it lacks
MISSING = {"zebrafish": {2, 3}, "mouse": {3}, "macaque": set(), "human": set()}
CENT = []
for k, (lab, key) in enumerate(SPECIES):
    draw_plane(Rk, k, key, 9 - k)
    cents = {}
    for j, (sname, sc, ca, cb) in enumerate(STATES):
        ca2, cb2 = ca + rs.normal(0, 0.025), cb + rs.normal(0, 0.04)
        cents[j] = Rk(k, ca2, cb2)
        if j in MISSING[key]:                                           # this species has no such state
            ax.add_patch(Ellipse(cents[j], 0.5, 0.2, fc="none", ec=sc, lw=1.2, ls=(0, (1.5, 1.5)), zorder=12))
            continue
        n = 28
        ua = np.clip(ca2 + rs.normal(0, 0.035, n), 0.03, 0.97)
        ub = np.clip(cb2 + rs.normal(0, 0.09, n), 0.06, 0.94)
        pts = np.array([Rk(k, u, v) for u, v in zip(ua, ub)])
        ax.scatter(pts[:, 0], pts[:, 1], s=8, color=sc, edgecolors="none", zorder=11)
    CENT.append(cents)
for ca, cb in zip(CENT, CENT[1:]):
    for j in range(len(STATES)):
        ax.plot([ca[j][0], cb[j][0]], [ca[j][1], cb[j][1]], color=MUTED, lw=0.9, ls=(0, (3, 2)), zorder=10)
header(ax, 6.85, 11.9, H3 - 0.18, "Mapped cell states")
fig.savefig(os.path.join(OUT, "transfer.png"), dpi=DPI, facecolor="white")
crop_to_ink(os.path.join(OUT, "transfer.png"))
plt.close(fig)
print("wrote", sorted(os.listdir(OUT)))
