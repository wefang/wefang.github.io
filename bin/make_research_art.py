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
# A small spatial branching process, in the manner of the research statement's simulations: cells
# live at continuous positions inside an organic tissue outline, divide after gamma-distributed
# waiting times, and scatter their daughters a little less each generation, so clones stay coherent
# while they intermix at their edges. At the end each cell reads its fate from a signal gradient
# plus a bias it inherited from its founder. The tissue is drawn as a Voronoi tiling of the cells.
from scipy.spatial import Voronoi
from matplotlib.patches import Polygon as Poly, Ellipse
import matplotlib.image as mpimg


def rad_of(th):
    return (1 + 0.10 * np.sin(2 * th + 0.5) + 0.07 * np.cos(3 * th) + 0.045 * np.sin(5 * th + 1.0)
            - 0.03 * np.cos(4 * th))


def inside(p):
    th = np.arctan2(p[1] - 0.5, p[0] - 0.5)
    return np.hypot((p[0] - 0.5) / 0.47, (p[1] - 0.5) / 0.45) <= 0.94 * rad_of(th)


def clip_in(p):
    c = np.array([0.5, 0.5])
    p = np.array(p, float)
    while not inside(p):
        p = c + 0.93 * (p - c)
    return p


rs = np.random.default_rng(21)
# a growing tissue in physical units (cell diameter 1): each generation every cell divides, the
# daughter is placed beside its mother, and overlapping cells push each other apart, so the tissue
# expands and each clone stays a coherent, ragged patch
GENS, CLONE_GEN, EARLY_GEN = 9, 3, 6
pos = np.zeros((1, 2))
clone = np.array([-1])
bias = np.zeros(1)
snap = {}


def relax(p, iters=30):
    for _ in range(iters):
        d = p[:, None, :] - p[None, :, :]
        r = np.hypot(d[..., 0], d[..., 1]) + np.eye(len(p)) * 9
        over = np.clip(1.0 - r, 0, None)
        p = p + 0.25 * (d / r[..., None] * over[..., None]).sum(1)
    return p


for g in range(1, GENS + 1):
    ang = rs.uniform(0, 2 * np.pi, len(pos))
    off = 0.5 * np.c_[np.cos(ang), np.sin(ang)]
    pos = np.vstack([pos - off, pos + off])
    clone = np.concatenate([clone, clone])
    bias = np.concatenate([bias, bias])
    if g == CLONE_GEN:                                   # label the founders of the clones we follow
        clone = np.arange(len(pos))
        bias = rs.normal(0, 0.09, len(pos))
    pos = relax(pos + rs.normal(0, 0.08, pos.shape))
    if g == EARLY_GEN:
        snap = dict(pos=pos.copy(), clone=clone.copy())
R_FINAL = np.max(np.hypot(*pos.T))


def to_tissue(p):
    """Map physical positions into the unit tissue frame, warped to the organic outline; the
    scale is fixed by the final tissue, so an earlier snapshot is drawn smaller."""
    th = np.arctan2(p[:, 1], p[:, 0])
    rr = np.hypot(*p.T) / R_FINAL
    warp = 0.94 * rad_of(th) * 0.97
    return np.c_[0.5 + 0.47 * warp * rr * np.cos(th), 0.5 + 0.45 * warp * rr * np.sin(th)]


P = to_tissue(pos)
mid = to_tissue(snap["pos"])
mid_clone = snap["clone"]
R_EARLY = np.max(np.hypot(*snap["pos"].T)) / R_FINAL
final = [dict(clone=int(c), fate=int(np.digitize(p[0] + b + rs.normal(0, 0.05), [0.40, 0.64])))
         for p, c, b in zip(P, clone, bias)]
alive_mid = [dict(clone=int(c)) for c in mid_clone]

CLONE_COLS = [TEAL, ORANGE, PURPLE, MAGENTA, BROWN, BLUE, "#66A61E", "#E6AB02"]
FATE_COLS = ["#F2C14E", "#5AA9A0", "#3D2C6B"]


def tissue(ax, x0, y0, w, h, pts, colors, title, scale=1.0):
    th = np.linspace(0, 2 * np.pi, 300)
    r = 0.94 * rad_of(th) * scale
    outline = np.c_[x0 + w * (0.5 + 0.47 * r * np.cos(th)), y0 + h * (0.5 + 0.45 * r * np.sin(th))]
    clip = Poly(outline, closed=True, fc="none", ec="none")
    ax.add_patch(clip)
    xy = np.c_[x0 + w * pts[:, 0], y0 + h * pts[:, 1]]
    far = np.array([[x0 - 5, y0 - 5], [x0 + w + 5, y0 - 5], [x0 - 5, y0 + h + 5], [x0 + w + 5, y0 + h + 5]])
    vor = Voronoi(np.vstack([xy, far]))
    for i in range(len(xy)):
        reg = vor.regions[vor.point_region[i]]
        if -1 in reg or not reg:
            continue
        tile = Poly(vor.vertices[reg], closed=True, fc=colors[i], ec="white", lw=0.9, zorder=3)
        ax.add_patch(tile)
        tile.set_clip_path(clip)
    ax.add_patch(Poly(outline, closed=True, fc="none", ec=MUTED, lw=1.1, zorder=4))
    text(ax, x0 + w / 2, y0 + h + 0.12, title, fontsize=12, color=INK)


fig, ax = canvas()
mid_cols = [CLONE_COLS[c["clone"] % 8] for c in alive_mid]
tissue(ax, 0.35, 0.55, 3.5, 2.7, mid, mid_cols, "Clones, early", scale=R_EARLY * 1.04)
ax.add_patch(FancyArrowPatch((3.85, 1.83), (4.3, 1.83), arrowstyle="-|>", mutation_scale=13, color=MUTED, lw=1.3))
tissue(ax, 4.45, 0.55, 3.5, 2.7, P, [CLONE_COLS[c["clone"] % 8] for c in final], "Clones, later")
tissue(ax, 8.25, 0.55, 3.5, 2.7, P, [FATE_COLS[c["fate"]] for c in final], "Fates")
grad = np.linspace(0, 1, 300).reshape(1, -1)
ax.imshow(grad, extent=(8.45, 11.55, 0.3, 0.44), aspect="auto",
          cmap=matplotlib.colors.LinearSegmentedColormap.from_list("sig", ["#FFFFFF", "#9E9AC8", "#3F007D"]))
text(ax, 8.45, 0.14, "Low signal", fontsize=10.5, ha="left")
text(ax, 11.55, 0.14, "High signal", fontsize=10.5, ha="right")
save(fig, "fate.png")

# =============================== 3. transfer ===============================
# Left half: a phylogeny with a silhouette per species, a functional genomics signal track for the
# same syntenic region in each, orthologous genes marked, and each organism's proposal onto human.
# Right half: the cell states of the same four species, aligned row to row. Human is sampled only
# to an earlier stage, so its later states are outlines.
SIL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "silhouettes")
W3, H3 = 12.0, 4.3
fig = plt.figure(figsize=(W3, H3), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W3)
ax.set_ylim(0, H3)
ax.set_aspect("equal")
ax.axis("off")
rs = np.random.default_rng(12)
SPECIES = [("Zebrafish", "zebrafish", BLUE), ("Mouse", "mouse", ORANGE), ("Macaque", "macaque", PURPLE), ("Human", "human", TEAL)]
ROW = [3.25, 2.4, 1.55, 0.7]                    # track baselines, human at the bottom


def hline(x0, x1, y):
    ax.plot([x0, x1], [y, y], color=INK, lw=1.5, solid_capstyle="round")


def vline(x, y0, y1):
    ax.plot([x, x], [y0, y1], color=INK, lw=1.5, solid_capstyle="round")


# phylogeny, root on the left: ((macaque, human), mouse), zebrafish
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
    silhouette(key, 1.35, y + 0.36, 0.58, 0.7, col)
    text(ax, 1.35, y - 0.1, lab, fontsize=10.5, color=col, fontweight="bold")

# signal tracks, generated from peaks placed per species; element positions drift between species
TX0, TX1 = 1.95, 6.15
xs = np.linspace(TX0, TX1, 600)
ELEM = {"zebrafish": [2.85, 3.85, 4.9], "mouse": [3.0, 4.05, 5.05], "macaque": [3.1, 4.15, 5.12], "human": [3.12, 4.17, 5.14]}
GA = {"zebrafish": 2.05, "mouse": 2.1, "macaque": 2.15, "human": 2.15}     # gene A, left flank
GB = {"zebrafish": 5.5, "mouse": 5.6, "macaque": 5.65, "human": 5.67}      # gene B, right flank


def gene(x, y, color, label=None):
    ax.add_patch(Polygon([(x, y - 0.07), (x + 0.26, y - 0.07), (x + 0.34, y), (x + 0.26, y + 0.07), (x, y + 0.07)],
                         closed=True, fc=color, ec="none", zorder=5))
    if label:
        text(ax, x + 0.17, y + 0.2, label, fontsize=9.5, color=INK, style="italic")


for (lab, key, col), y in zip(SPECIES, ROW):
    sig = 0.03 * rs.random(xs.size)
    for k, e in enumerate(ELEM[key]):
        h = 0.52 if not (key == "human" and k == 0) else 0.08          # element 1 is silent in human
        sig += h * np.exp(-((xs - e) ** 2) / (2 * 0.05 ** 2))
    for e in rs.uniform(TX0 + 0.5, TX1 - 0.7, 3):                        # species-specific peaks
        sig += 0.16 * np.exp(-((xs - e) ** 2) / (2 * 0.035 ** 2))
    ax.fill_between(xs, y, y + sig, color=col, alpha=0.9, lw=0, zorder=3)
    ax.plot([TX0, TX1], [y, y], color=GREY, lw=1, zorder=2)
    top = key == "zebrafish"
    gene(GA[key], y, INK, "Gene A" if top else None)
    gene(GB[key], y, INK, "Gene B" if top else None)
# orthologous genes joined across species, syntenic regions shaded between neighbouring tracks
for (_, a, _), (_, b, _), ya, yb in zip(SPECIES, SPECIES[1:], ROW, ROW[1:]):
    for g in (GA, GB):
        ax.plot([g[a] + 0.17, g[b] + 0.17], [ya - 0.08, yb + 0.08], color=MUTED, lw=1.0, ls=(0, (2, 2)), zorder=2)
    for k in range(3):
        xa, xb = ELEM[a][k], ELEM[b][k]
        ax.add_patch(Polygon([(xa - 0.17, ya - 0.01), (xa + 0.17, ya - 0.01), (xb + 0.17, yb + 0.54), (xb - 0.17, yb + 0.54)],
                             closed=True, fc="#DAD5E0", ec="none", alpha=0.85, zorder=1))
text(ax, (GA["human"] + GB["human"]) / 2 + 0.17, ROW[3] - 0.25, "Orthologous genes joined, syntenic regions shaded",
     fontsize=9.5, color=MUTED)
# each model organism's proposal for element 2, onto the human track
for (lab, key, col), y in zip(SPECIES[:3], ROW[:3]):
    ax.add_patch(FancyArrowPatch((ELEM[key][1] + 0.12, y + 0.28), (ELEM["human"][1] + 0.1, ROW[3] + 0.5),
                                 connectionstyle="arc3,rad=-0.3", arrowstyle="-|>", mutation_scale=11,
                                 color=col, lw=1.4, zorder=6))
text(ax, (TX0 + TX1) / 2, 4.08, "Regulatory elements", fontsize=12, color=INK)

# right half: cell states, one small branching trajectory per species, aligned across species
CX0 = 6.75
CT = {"prog": "#8C8C8C", "a": "#E6AB02", "b": "#E7298A"}


def mix(c0, c1, f):
    a0, a1 = np.array(matplotlib.colors.to_rgb(c0)), np.array(matplotlib.colors.to_rgb(c1))
    return a0 + (a1 - a0) * f


def branch_pts(ox, oy, human):
    """A cloud of cells along a stem that forks in two, colored from progenitor into each fate.
    Returns anchor points for the alignment lines."""
    anchors = {}
    n = 70
    t = np.sort(rs.random(n))
    pts = np.c_[ox + 1.6 * t, oy + rs.normal(0, 0.05, n)]
    ax.scatter(pts[:, 0], pts[:, 1], s=9, color=CT["prog"], edgecolors="none", alpha=0.9, zorder=4)
    anchors["fork"] = np.array([ox + 1.6, oy])
    for tag, sgn in (("a", 1), ("b", -1)):
        u = np.sort(rs.random(90))
        if human:
            u = u[u < 0.42]                                         # later states not yet sampled
        spread = 0.04 + 0.03 * u
        px = ox + 1.6 + 2.9 * u
        py = oy + sgn * 0.32 * (3 * u ** 2 - 2 * u ** 3) + rs.normal(0, 1, u.size) * spread
        cols = [mix(CT["prog"], CT[tag], min(1, 0.25 + 1.3 * v)) for v in u]
        ax.scatter(px, py, s=9, c=cols, edgecolors="none", alpha=0.9, zorder=4)
        end = np.array([ox + 4.5, oy + sgn * 0.32])
        if human:
            uu = np.linspace(0.42, 1, 30)
            ax.plot(ox + 1.6 + 2.9 * uu, oy + sgn * 0.32 * (3 * uu ** 2 - 2 * uu ** 3), color=GREY, lw=1.1,
                    ls=(0, (2, 2)), zorder=2)
            ax.add_patch(Ellipse(end, 0.6, 0.22, fc="none", ec=CT[tag], lw=1.2, ls=(0, (2, 1.5)), zorder=3))
        anchors[tag] = end
        anchors[tag + "_mid"] = np.array([ox + 1.6 + 2.9 * 0.4, oy + sgn * 0.32 * (3 * 0.16 - 2 * 0.064)])
    return anchors


ANCH = [branch_pts(CX0, y + 0.3, key == "human") for (lab, key, col), y in zip(SPECIES, ROW)]
for ta, tb in zip(ANCH, ANCH[1:]):
    for k in ("fork", "a_mid", "b_mid", "a", "b"):
        ax.plot([ta[k][0], tb[k][0]], [ta[k][1], tb[k][1]], color="#7FA6D6", lw=0.9, ls=(0, (3, 2)), zorder=1)
text(ax, CX0 + 2.25, 4.08, "Cell states, aligned across species", fontsize=12, color=INK)
for k, (tag, lab) in enumerate((("prog", "Progenitor"), ("a", "Fate A"), ("b", "Fate B"))):
    ax.scatter([CX0 + 0.1 + k * 1.55], [0.12], s=30, color=CT[tag], edgecolors="none")
    text(ax, CX0 + 0.22 + k * 1.55, 0.12, lab, fontsize=10, ha="left")
fig.savefig(os.path.join(OUT, "transfer.png"), dpi=DPI, facecolor="white")
plt.close(fig)
print("wrote", sorted(os.listdir(OUT)))
