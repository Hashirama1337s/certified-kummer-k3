"""Draw the README figures of certified-kummer-k3 from the data the repository ships.

Usage:  python figures/make_figures.py <repo root> [--png]
Needs Python 3.11 and matplotlib 3.11 (nothing else). Writes figures/*.svg next to this script
(or into ./figures when the script does not sit in a folder named figures). --png also writes PNG previews.

Every plotted number is read from shipped files:
  results/wp5/result_wp5.json   stageB: for each of the 30 bubble-size intervals [a_lo, a_hi] the certified upper
                                bounds d2 (closeness) and smallness U (continuity step), and a0      [sealed test WP5]
  code/wp5/wp5.py               the closing rule of the continuity step (smallness < 2), checked as text below
"""
import json
import os
import re
import sys
from fractions import Fraction as F

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

# ---------------------------------------------------------------- brand style (figures/STYLE.md)
BG, BORDER, GRID, AXC, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]
W_IN = 12.5  # 900 pt viewBox width

plt.rcParams.update({
    "svg.fonttype": "none", "svg.hashsalt": "moki-julio",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 16, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "axes.facecolor": BG, "figure.facecolor": "none", "savefig.facecolor": "none",
    "axes.edgecolor": AXC, "axes.labelcolor": TXT, "xtick.color": AXC, "ytick.color": AXC,
    "text.color": TXT, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.fontsize": 14, "lines.linewidth": 2.0, "axes.axisbelow": True,
})
PNG = "--png" in sys.argv


def new_fig(h_in):
    fig = plt.figure(figsize=(W_IN, h_in))
    fig.patches.append(FancyBboxPatch((0.004, 0.006), 0.992, 0.988, boxstyle="round,pad=0,rounding_size=0.015",
                                      transform=fig.transFigure, facecolor=BG, edgecolor=BORDER, linewidth=1.5,
                                      zorder=-10))
    return fig


def title(fig, main, sub):
    h = fig.get_figheight()
    fig.text(0.035, 1 - 0.42 / h, main, fontsize=21, fontweight="bold", color=TXT, va="top")
    fig.text(0.035, 1 - 0.88 / h, sub, fontsize=15, color=AXC, va="top")


def save(fig, out, name):
    path = os.path.join(out, name)
    fig.savefig(path, format="svg", metadata={"Date": None, "Creator": None})
    with open(path, encoding="utf-8") as fh:  # mathtext pieces name one font only; give them the full stack
        svg = fh.read().replace("font-family: 'Segoe UI';", "font-family: 'Segoe UI', 'Helvetica', 'Arial', sans-serif;")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(svg)
    if PNG:
        prev = os.path.join(out, os.pardir, "_preview")
        os.makedirs(prev, exist_ok=True)
        fig.savefig(os.path.join(prev, name.replace(".svg", ".png")), format="png", dpi=144)
    plt.close(fig)
    print("wrote figures/" + name)


def up(x, nd):
    """float -> decimal string rounded UP at nd decimals (exact, via the float's binary value)."""
    q = F(x) * 10 ** nd
    n = -((-q.numerator) // q.denominator)
    s = str(n).rjust(nd + 1, "0")
    return s[:-nd] + "." + s[-nd:]


def stage_b(root):
    with open(os.path.join(root, "results/wp5/result_wp5.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    rows = []
    for key, v in d["stageB"].items():
        lo, hi = (float(t) for t in re.findall(r"[0-9.]+", key))
        assert v["closes"] in (True, "True")
        rows.append((lo, hi, float(v["d2"]), float(v["smallness"])))
    rows.sort()
    assert len(rows) == 30 and rows[0][0] == 0 and abs(rows[-1][1] - float(d["a0"])) < 1e-15
    assert all(abs(a[1] - b[0]) < 1e-15 for a, b in zip(rows, rows[1:]))  # the intervals tile (0, a0]
    return rows, float(d["a0"])


# ---------------------------------------------------------------- figure 1: the closeness bound (hero)
def fig_bound(root, out):
    rows, a0 = stage_b(root)
    fig = new_fig(7.4)
    title(fig, "How close is the hand-built K3 geometry?",
          "Certified upper bound d on how far the perfect (Ricci-flat) geometry can differ from the glued one")
    ax = fig.add_axes([0.10, 0.13, 0.86, 0.66])
    floor = 5e-7
    for lo, hi, d2, _ in rows:
        ax.add_patch(Rectangle((lo, floor), hi - lo, d2 - floor, facecolor=CAT[0], alpha=0.28, lw=0))
        ax.plot([lo, hi], [d2, d2], color=CAT[0], lw=2.2, solid_capstyle="butt")
    for a, b in zip(rows, rows[1:]):
        ax.plot([a[1], a[1]], [a[2], b[2]], color=CAT[0], lw=1)
    ax.set_yscale("log")
    ax.set_xlim(0, a0 * 1.13)
    ax.set_ylim(floor, 0.6)
    ax.axvspan(a0, a0 * 1.13, color=AXC, alpha=0.10, lw=0)
    ax.text(a0 * 1.065, 1.5e-6, "not\ncovered\n(a > %s)" % a0, color=AXC, fontsize=13, ha="center")
    ax.set_xlabel("bubble size a (size of the 16 smoothing bubbles)")
    ax.set_ylabel("closeness bound d (log scale)")
    picks = {}
    for lo, hi, d2, _ in rows:
        if abs(hi - 0.01) < 1e-12:
            picks["0.01"] = (hi, d2)
        if abs(hi - a0) < 1e-12:
            picks["a0"] = (hi, d2)
    h, d = picks["0.01"]
    ax.plot(h, d, "o", ms=9, color=ACC, zorder=5)
    ax.annotate("a = 0.01:  d ≤ %s\n(under 0.3 percent)" % up(d, 6), (h, d), xytext=(0.0115, 2.2e-4),
                color=TXT, fontsize=14, va="center", arrowprops=dict(arrowstyle="-", color=AXC, lw=1.2))
    h, d = picks["a0"]
    ax.plot(h, d, "o", ms=9, color=ACC, zorder=5)
    ax.annotate("a = %s:  d ≤ %s" % (a0, up(d, 6)), (h, d), xytext=(0.0175, 0.25), color=TXT, fontsize=14,
                va="center", arrowprops=dict(arrowstyle="-", color=AXC, lw=1.2))
    ax.text(0.0012, 3e-2, "each step is a proved upper bound for every a under it:\nthe true d(a) lies on or below the step (shaded)",
            color=AXC, fontsize=14, va="center")
    save(fig, out, "closeness-bound.svg")
    return dict(n=len(rows), first=rows[0], at001=picks["0.01"], at_a0=picks["a0"],
                printed=(up(picks["0.01"][1], 6), up(picks["a0"][1], 6)))


# ---------------------------------------------------------------- figure 2: room to spare at a = 0.03
def threshold(root):
    src = open(os.path.join(root, "code/wp5/wp5.py"), encoding="utf-8").read()
    m = re.search(r"r\['smallness'\] < (\d+)", src)
    assert m, "closing rule not found in code/wp5/wp5.py"
    return int(m.group(1))


def fig_smallness(root, out):
    rows, a0 = stage_b(root)
    lim = threshold(root)
    fig = new_fig(6.6)
    title(fig, "Where the certificate stops, and why",
          "The proof's safety number U must stay below %d; at the last bubble size it is still far below" % lim)
    ax = fig.add_axes([0.10, 0.14, 0.86, 0.62])
    for lo, hi, _, u in rows:
        ax.add_patch(Rectangle((lo, 0), hi - lo, u, facecolor=ACC, alpha=0.22, lw=0))
        ax.plot([lo, hi], [u, u], color=ACC, lw=2.2, solid_capstyle="butt")
    ax.axhline(lim, color=CAT[1], lw=2)
    ax.text(0.0005, lim + 0.04, "the method needs U < %d" % lim, color=CAT[1], fontsize=14, va="bottom")
    lo, hi, _, u = rows[-1]
    ax.annotate("at a = %s:  U ≤ %s" % (a0, up(u, 4)), (hi, u), xytext=(0.0185, 1.55), color=TXT,
                fontsize=14, va="center", arrowprops=dict(arrowstyle="-", color=AXC, lw=1.2))
    ax.text(0.0005, 0.55, "U for each bubble-size interval\n(certified upper bound)", color=AXC, fontsize=14)
    ax.set_xlim(0, a0 * 1.3)
    ax.set_ylim(0, lim * 1.15)
    ax.axvspan(a0, a0 * 1.3, color=AXC, alpha=0.10, lw=0)
    ax.text(a0 * 1.15, 0.62, "stops at a = %s:\na lemma about the\nsurface's diameter\nis proved only\nthat far" % a0,
            color=TXT, fontsize=13, ha="center", va="center")
    ax.set_xlabel("bubble size a")
    ax.set_ylabel("safety number U")
    save(fig, out, "where-it-stops.svg")
    return dict(threshold=lim, U_first=rows[0][3], U_last=rows[-1][3], printed=up(rows[-1][3], 4))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        sys.exit(__doc__)
    root = args[0]
    here = os.path.dirname(os.path.abspath(__file__))
    out = here if os.path.basename(here) == "figures" else os.path.join(here, "figures")
    os.makedirs(out, exist_ok=True)
    info = dict(bound=fig_bound(root, out), smallness=fig_smallness(root, out))
    if "--info" in sys.argv:
        print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
