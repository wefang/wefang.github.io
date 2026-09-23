"""Draw the three research-area images for the about page, wide and short so each spans the page.

  lineage.png   a dated lineage tree with the recorder edits that happened on its branches, and the
                character matrix the sampled cells carry (colored as in the research statement)
  fate.png      a small simulation: clones grow and divide across a tissue with a signal gradient,
                shown once colored by clone and once colored by the fate each cell took
  transfer.png  a phylogeny joining model organisms and human, a functional genomics signal track
                for each, syntenic blocks between them, and each organism's proposal onto human

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
fig, ax = canvas()
rng = np.random.default_rng(4)
NSITE, NTIP = 8, 12
# a random binary tree with dated nodes, drawn sideways: time runs left to right
X0, X1 = 0.55, 5.55                      # root time, sampling time
Y0, Y1 = 0.45, 3.0                      # tip rows span


def split(tips, t, depth):
    """Divide a set of tips at time t; returns nested (time, children) structure."""
    if len(tips) == 1:
        return ("tip", tips[0])
    k = int(rng.integers(1, len(tips)))
    k = max(1, min(len(tips) - 1, len(tips) // 2 + int(rng.integers(-1, 2))))
    tl = t + (X1 - t) * rng.uniform(0.18, 0.42)
    tr = t + (X1 - t) * rng.uniform(0.18, 0.42)
    return ("node", t, split(tips[:k], tl, depth + 1), split(tips[k:], tr, depth + 1))


tree = split(list(range(NTIP)), X0, 0)
tip_y = {i: Y1 - i * (Y1 - Y0) / (NTIP - 1) for i in range(NTIP)}
states = {i: [None] * NSITE for i in range(NTIP)}
used_sites = []


def ypos(n):
    if n[0] == "tip":
        return tip_y[n[1]]
    return (ypos(n[2]) + ypos(n[3])) / 2


def tips_of(n):
    return [n[1]] if n[0] == "tip" else tips_of(n[2]) + tips_of(n[3])


def draw(n, parent_t, record):
    """Draw the branch into n from parent_t; an edit may land on it and is inherited below."""
    y = ypos(n)
    t = X1 if n[0] == "tip" else n[1]
    ax.plot([parent_t, t], [y, y], color=INK, lw=1.6, solid_capstyle="round", zorder=2)
    record = list(record)
    free = [s for s in range(NSITE) if record[s] is None]
    if parent_t < t - 0.35 and free and rng.random() < 0.8:
        s = int(rng.choice(free))
        a = ALLELES[int(rng.integers(len(ALLELES)))]
        record[s] = a
        ex = parent_t + (t - parent_t) * rng.uniform(0.35, 0.65)
        ax.add_patch(Rectangle((ex - 0.08, y - 0.08), 0.16, 0.16, fc=a, ec="white", lw=0.8, zorder=4))
    if n[0] == "tip":
        states[n[1]] = record
        return
    ya, yb = ypos(n[2]), ypos(n[3])
    ax.plot([t, t], [ya, yb], color=INK, lw=1.6, solid_capstyle="round", zorder=2)
    draw(n[2], t, record)
    draw(n[3], t, record)


ax.plot([X0 - 0.35, X0], [ypos(tree)] * 2, color=INK, lw=1.6, zorder=2)
ax.plot([X0], [ypos(tree)], "o", ms=5, color=INK, zorder=3)
draw(tree[2], X0, [None] * NSITE)
draw(tree[3], X0, [None] * NSITE)
ax.plot([X0, X0], [ypos(tree[2]), ypos(tree[3])], color=INK, lw=1.6, zorder=2)
for i in range(NTIP):
    ax.plot([X1], [tip_y[i]], "o", ms=4.2, color=INK, zorder=3)
# time axis
ax.add_patch(FancyArrowPatch((X0 - 0.3, 0.14), (X1, 0.14), arrowstyle="-|>", mutation_scale=11, color=GREY, lw=1.1))
text(ax, X1 + 0.12, 0.14, "Time", fontsize=11, ha="left")

# the character matrix the sampled cells carry: rows are cells, columns are recorder sites
MX = 6.35
cw, ch = 0.42, (Y1 - Y0) / (NTIP - 1) * 0.84
for i in range(NTIP):
    rec = states[i]
    for s in range(NSITE):
        c = rec[s] if rec[s] is not None else INK
        if rng.random() < 0.07:
            c = DROPOUT                                    # the readout failed at this site
        ax.add_patch(Rectangle((MX + s * (cw + 0.05), tip_y[i] - ch / 2), cw, ch, fc=c, ec="none", zorder=3))
    ax.plot([X1 + 0.1, MX - 0.08], [tip_y[i]] * 2, color=HAIR, lw=0.9, ls=(0, (2, 2)), zorder=1)
text(ax, MX + NSITE * (cw + 0.05) / 2, 3.33, "Recorded sites", fontsize=12, color=INK)
text(ax, (X0 + X1) / 2, 3.33, "Lineage with recording events", fontsize=12, color=INK)
# legend
LX = MX + NSITE * (cw + 0.05) + 0.25
for k, (c, lab) in enumerate(((INK, "Unedited"), (DROPOUT, "Dropout"), (TEAL, "Edits"))):
    y = 2.75 - k * 0.42
    if lab == "Edits":
        for j, cc in enumerate(ALLELES[:4]):
            ax.add_patch(Rectangle((LX + j * 0.1, y - 0.1), 0.09, 0.2, fc=cc, ec="none"))
    else:
        ax.add_patch(Rectangle((LX, y - 0.1), 0.37, 0.2, fc=c, ec="none"))
    text(ax, LX + 0.5, y, lab, fontsize=10.5, ha="left")
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
BLOCKS = {  # (start, end, color, strand) along the genome axis, in the order they occur
    "zebrafish": [(2.45, 3.35, "#C9B8E6", 1), (3.45, 4.75, "#F6C9A0", -1), (4.8, 5.35, "#B8D8C8", 1)],
    "mouse": [(2.45, 3.4, "#C9B8E6", 1), (3.5, 4.55, "#F6C9A0", 1), (4.65, 5.45, "#B8D8C8", 1)],
    "macaque": [(2.5, 3.45, "#C9B8E6", 1), (3.55, 4.55, "#F6C9A0", 1), (4.62, 5.4, "#B8D8C8", 1)],
    "human": [(2.5, 3.45, "#C9B8E6", 1), (3.55, 4.57, "#F6C9A0", 1), (4.63, 5.42, "#B8D8C8", 1)],
}


def gene(x, y, color, label=None):
    ax.add_patch(Polygon([(x, y - 0.07), (x + 0.26, y - 0.07), (x + 0.34, y), (x + 0.26, y + 0.07), (x, y + 0.07)],
                         closed=True, fc=color, ec="none", zorder=5))
    if label:
        text(ax, x + 0.17, y + 0.2, label, fontsize=9.5, color=INK, style="italic")


for (lab, key, col), y in zip(SPECIES, ROW):
    # the synteny blocks, drawn as a band just under the track with their orientation
    for s0, s1, bc, strand in BLOCKS[key]:
        ax.add_patch(Rectangle((s0, y - 0.2), s1 - s0, 0.1, fc=bc, ec="none", zorder=2))
        ax.add_patch(FancyArrowPatch((s0 + 0.08, y - 0.15) if strand > 0 else (s1 - 0.08, y - 0.15),
                                     (s1 - 0.08, y - 0.15) if strand > 0 else (s0 + 0.08, y - 0.15),
                                     arrowstyle="-|>", mutation_scale=6, color=MUTED, lw=0.7, zorder=3))
    sig = 0.03 * rs.random(xs.size)
    for k, e in enumerate(ELEM[key]):
        h = 0.5 if not (key == "human" and k == 0) else 0.07           # element 1 is silent in human
        sig += h * np.exp(-((xs - e) ** 2) / (2 * 0.045 ** 2))
    for e in rs.uniform(TX0 + 0.5, TX1 - 0.8, 3):                       # species-specific peaks
        sig += 0.15 * np.exp(-((xs - e) ** 2) / (2 * 0.03 ** 2))
    ax.fill_between(xs, y, y + sig, color=col, alpha=0.9, lw=0, zorder=3)
    ax.plot([TX0, TX1], [y, y], color=GREY, lw=1, zorder=2)
    top = key == "zebrafish"
    gene(GA[key], y, INK, "Gene A" if top else None)
    gene(GB[key], y, INK, "Gene B" if top else None)
# orthologous genes joined across species
for (_, a, _), (_, b, _), ya, yb in zip(SPECIES, SPECIES[1:], ROW, ROW[1:]):
    for g in (GA, GB):
        ax.plot([g[a] + 0.17, g[b] + 0.17], [ya - 0.22, yb + 0.1], color=MUTED, lw=0.9, ls=(0, (2, 2)), zorder=1)
# each model organism's proposal for the element in the second block, onto the human track
for (lab, key, col), y in zip(SPECIES[:3], ROW[:3]):
    ax.add_patch(FancyArrowPatch((ELEM[key][1] + 0.1, y + 0.3), (ELEM["human"][1] + 0.08, ROW[3] + 0.5),
                                 connectionstyle="arc3,rad=-0.3", arrowstyle="-|>", mutation_scale=11,
                                 color=col, lw=1.4, zorder=6))
text(ax, (TX0 + TX1) / 2, 4.38, "Regulatory elements", fontsize=12, color=INK)
text(ax, 4.0, 0.2, "Synteny blocks shown under each track, with orientation", fontsize=9.5, color=MUTED)

# ---------------- right: a stack of embedding planes, one per species ----------------
STAGE = plt.cm.YlGnBu


def Pk(k, a, b):
    ox, oy = 6.55 + 0.34 * k, ROW[k] - 0.08
    return ox + 3.5 * a + 0.7 * b, oy + 0.62 * b


def arc(k, p0, p1, t0, t1, upto=1.0, bow=0.12):
    mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    ctrl = (mx - bow * dy, my + bow * dx)
    us = np.linspace(0, 1, 26)
    for u in us:
        a = (1 - u) ** 2 * p0[0] + 2 * (1 - u) * u * ctrl[0] + u ** 2 * p1[0]
        b = (1 - u) ** 2 * p0[1] + 2 * (1 - u) * u * ctrl[1] + u ** 2 * p1[1]
        if u > upto:
            continue
        for _ in range(3):
            aa = min(max(a + rs.normal(0, 0.018), 0.03), 0.97)
            bb = min(max(b + rs.normal(0, 0.05), 0.05), 0.95)
            x, y = Pk(k, aa, bb)
            ax.plot([x], [y], ".", ms=3.0, color=STAGE(0.2 + 0.75 * (t0 + u * (t1 - t0))), zorder=10 - k)
    return Pk(k, *p1)


TIPS = []
for k, (lab, key, col) in enumerate(SPECIES):
    ax.add_patch(Polygon([Pk(k, 0, 0), Pk(k, 1, 0), Pk(k, 1, 1), Pk(k, 0, 1)], closed=True,
                         fc="#FAFAFB", ec=col, lw=1.1, alpha=0.95, zorder=9 - k))
    human = key == "human"
    fork = arc(k, (0.06, 0.45), (0.42, 0.5), 0.0, 0.42)
    up = arc(k, (0.42, 0.5), (0.95, 0.85), 0.42, 1.0, upto=0.45 if human else 1.0)
    dn = arc(k, (0.42, 0.5), (0.95, 0.18), 0.42, 1.0, upto=0.45 if human else 1.0, bow=-0.12)
    if human:
        for p in (up, dn):
            ax.add_patch(Ellipse(p, 0.42, 0.16, fc="none", ec=col, lw=1.1, ls=(0, (1.5, 1.5)), zorder=11))
    TIPS.append((fork, up, dn))
# corresponding states linked from plane to plane
for (fa, ua, da), (fb, ub, db) in zip(TIPS, TIPS[1:]):
    for p, q in ((fa, fb), (ua, ub), (da, db)):
        ax.plot([p[0], q[0]], [p[1], q[1]], color="#7FA6D6", lw=0.9, ls=(0, (3, 2)), zorder=1)
text(ax, 8.9, 4.38, "Cell states, aligned across species", fontsize=12, color=INK)
for i in range(6):
    ax.add_patch(Rectangle((7.2 + i * 0.2, 0.14), 0.19, 0.11, fc=STAGE(0.2 + 0.75 * i / 5), ec="none"))
text(ax, 7.1, 0.2, "Stage", fontsize=10, ha="right")
text(ax, 8.55, 0.2, "early to late", fontsize=9.5, ha="left")
fig.savefig(os.path.join(OUT, "transfer.png"), dpi=DPI, facecolor="white")
plt.close(fig)
print("wrote", sorted(os.listdir(OUT)))
