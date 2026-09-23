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
text(ax, X1 + 0.12, 0.14, "time", fontsize=11, ha="left")

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
text(ax, MX + NSITE * (cw + 0.05) / 2, 3.33, "recorded sites", fontsize=12, color=INK)
text(ax, (X0 + X1) / 2, 3.33, "lineage with recording events", fontsize=12, color=INK)
# legend
LX = MX + NSITE * (cw + 0.05) + 0.25
for k, (c, lab) in enumerate(((INK, "unedited"), (DROPOUT, "dropout"), (TEAL, "edits"))):
    y = 2.75 - k * 0.42
    if lab == "edits":
        for j, cc in enumerate(ALLELES[:4]):
            ax.add_patch(Rectangle((LX + j * 0.1, y - 0.1), 0.09, 0.2, fc=cc, ec="none"))
    else:
        ax.add_patch(Rectangle((LX, y - 0.1), 0.37, 0.2, fc=c, ec="none"))
    text(ax, LX + 0.5, y, lab, fontsize=10.5, ha="left")
save(fig, "lineage.png")

# =============================== 2. fate ===============================
# A small agent simulation on a grid: founders seed a tissue, clones grow by division into
# free neighbouring sites, and each cell's fate is read off the local signal plus a little
# inherited bias. Then the same tissue is drawn twice, by clone and by fate.
rng = np.random.default_rng(7)
GX, GY = 46, 22
clone = -np.ones((GY, GX), dtype=int)
bias = np.zeros((GY, GX))
founders = [(int(rng.integers(2, GY - 2)), int(x)) for x in np.linspace(3, GX - 4, 9)]
frontier = []
for k, (y, x) in enumerate(founders):
    clone[y, x] = k
    bias[y, x] = rng.normal(0, 0.10)
    frontier.append((y, x))
nbrs = [(-1, 0), (1, 0), (0, -1), (0, 1)]
while frontier:
    i = int(rng.integers(len(frontier)))
    y, x = frontier[i]
    free = [(y + dy, x + dx) for dy, dx in nbrs
            if 0 <= y + dy < GY and 0 <= x + dx < GX and clone[y + dy, x + dx] < 0]
    if not free:
        frontier.pop(i)
        continue
    ny, nx = free[int(rng.integers(len(free)))]
    clone[ny, nx] = clone[y, x]
    bias[ny, nx] = bias[y, x] + rng.normal(0, 0.03)       # a daughter inherits its mother's bias
    frontier.append((ny, nx))
signal = np.tile(np.linspace(0, 1, GX), (GY, 1))
score = signal + bias + rng.normal(0, 0.06, size=signal.shape)
fate = np.digitize(score, [0.36, 0.68])                      # three fates along the gradient

CLONE_COLS = [TEAL, ORANGE, PURPLE, MAGENTA, BROWN, BLUE, "#66A61E", "#E6AB02", "#666666"]
FATE_COLS = ["#F3D9A4", "#E7298A", "#4B2E83"]
fig, ax = canvas()
cs = 0.112                                                 # cell size
panels = ((0.3, "clones", lambda y, x: CLONE_COLS[clone[y, x] % len(CLONE_COLS)]),
          (6.25, "fates", lambda y, x: FATE_COLS[fate[y, x]]))
for ox, title, colf in panels:
    oy = 0.75
    for y in range(GY):
        for x in range(GX):
            ax.add_patch(Rectangle((ox + x * cs, oy + y * cs), cs * 0.92, cs * 0.92, fc=colf(y, x), ec="none"))
    # clone outlines on both panels, so the fate panel still shows who is related to whom
    for y in range(GY):
        for x in range(GX):
            if x + 1 < GX and clone[y, x] != clone[y, x + 1]:
                ax.plot([ox + (x + 1) * cs - 0.005] * 2, [oy + y * cs, oy + (y + 1) * cs], color="white", lw=1.4)
            if y + 1 < GY and clone[y, x] != clone[y + 1, x]:
                ax.plot([ox + x * cs, ox + (x + 1) * cs], [oy + (y + 1) * cs - 0.005] * 2, color="white", lw=1.4)
    text(ax, ox + GX * cs / 2, oy + GY * cs + 0.25, title, fontsize=13, color=INK)
    grad = np.linspace(0, 1, 300).reshape(1, -1)
    ax.imshow(grad, extent=(ox, ox + GX * cs, 0.32, 0.52), aspect="auto",
              cmap=matplotlib.colors.LinearSegmentedColormap.from_list("sig", ["#FFFFFF", "#8C6BB1", "#3F007D"]))
    text(ax, ox - 0.05, 0.42, "", fontsize=1)
text(ax, 0.3, 0.12, "low signal", fontsize=10.5, ha="left")
text(ax, 0.3 + GX * cs, 0.12, "high signal", fontsize=10.5, ha="right")
text(ax, 6.25, 0.12, "low signal", fontsize=10.5, ha="left")
text(ax, 6.25 + GX * cs, 0.12, "high signal", fontsize=10.5, ha="right")
save(fig, "fate.png")

# =============================== 3. transfer ===============================
fig, ax = canvas()
rng = np.random.default_rng(12)
SPECIES = [("zebrafish", BLUE), ("mouse", ORANGE), ("macaque", PURPLE), ("human", TEAL)]
TY = [2.72, 1.92, 1.12, 0.32]                   # track baselines, human at the bottom
TX0, TX1 = 3.05, 11.7
# the phylogeny, sideways, root on the left
PX = 0.25
ax.plot([PX, 0.7], [1.5, 1.5], color=INK, lw=1.6)
ax.plot([0.7, 0.7], [TY[0] + 0.25, 0.95], color=INK, lw=1.6)
ax.plot([0.7, 1.55], [TY[0] + 0.25, TY[0] + 0.25], color=INK, lw=1.6)      # zebrafish
ax.plot([0.7, 1.05], [0.95, 0.95], color=INK, lw=1.6)
ax.plot([1.05, 1.05], [TY[1] + 0.25, 0.6], color=INK, lw=1.6)
ax.plot([1.05, 1.55], [TY[1] + 0.25, TY[1] + 0.25], color=INK, lw=1.6)      # mouse
ax.plot([1.05, 1.3], [0.6, 0.6], color=INK, lw=1.6)
ax.plot([1.3, 1.3], [TY[2] + 0.25, TY[3] + 0.25], color=INK, lw=1.6)
ax.plot([1.3, 1.55], [TY[2] + 0.25, TY[2] + 0.25], color=INK, lw=1.6)      # macaque
ax.plot([1.3, 1.55], [TY[3] + 0.25, TY[3] + 0.25], color=INK, lw=1.6)      # human
for (name, col), y in zip(SPECIES, TY):
    text(ax, 1.65, y + 0.25, name, fontsize=12, color=col, ha="left", fontweight="bold")

# a regulatory element and its region, placed per species; the positions drift between species
# (synteny holds, spacing does not), and one element is repurposed in human
ELEM = {"zebrafish": [4.4, 6.6, 9.3], "mouse": [4.9, 7.3, 9.9], "macaque": [5.1, 7.5, 10.1], "human": [5.2, 7.6, 10.2]}
GENE = {"zebrafish": 10.6, "mouse": 11.0, "macaque": 11.05, "human": 11.1}
xs = np.linspace(TX0, TX1, 700)


def track(y, name, col):
    base = 0.04 * rng.random(xs.size)
    sig = base.copy()
    for k, e in enumerate(ELEM[name]):
        h = 0.5 if not (name == "human" and k == 0) else 0.1      # element 1 is quiet in human
        sig += h * np.exp(-((xs - e) ** 2) / (2 * 0.09 ** 2))
    for e in rng.uniform(TX0, TX1, 4):                             # a few species-specific peaks
        sig += 0.15 * np.exp(-((xs - e) ** 2) / (2 * 0.06 ** 2))
    ax.fill_between(xs, y, y + sig, color=col, alpha=0.85, lw=0)
    ax.plot([TX0, TX1], [y, y], color=HAIR, lw=1)
    g = GENE[name]
    ax.add_patch(Polygon([(g, y - 0.09), (g + 0.42, y - 0.09), (g + 0.52, y), (g + 0.42, y + 0.09), (g, y + 0.09)],
                         closed=True, fc=INK, ec="none", zorder=3))


for (name, col), y in zip(SPECIES, TY):
    track(y, name, col)
# syntenic blocks between neighbouring tracks, shaded bands joining the same region
for (a, ya), (b, yb) in zip(zip([s[0] for s in SPECIES], TY), zip([s[0] for s in SPECIES][1:], TY[1:])):
    for k in range(3):
        xa, xb = ELEM[a][k], ELEM[b][k]
        ax.add_patch(Polygon([(xa - 0.3, ya - 0.01), (xa + 0.3, ya - 0.01), (xb + 0.3, yb + 0.56), (xb - 0.3, yb + 0.56)],
                             closed=True, fc="#D9D3DE", ec="none", alpha=0.85, zorder=0))
# each organism's proposal for element 2, drawn onto the human track
for (name, col), y in zip(SPECIES[:3], TY[:3]):
    x0 = ELEM[name][1] + 0.35
    ax.add_patch(FancyArrowPatch((x0, y + 0.25), (ELEM["human"][1] + 0.12, TY[3] + 0.52),
                                 connectionstyle="arc3,rad=-0.32", arrowstyle="-|>", mutation_scale=12,
                                 color=col, lw=1.5, zorder=4))
text(ax, ELEM["human"][1] - 0.35, 0.78, "proposed human element", fontsize=11, color=INK, ha="right")
save(fig, "transfer.png")
print("wrote", sorted(os.listdir(OUT)))
