#!/usr/bin/env python3
"""Figures and markdown tables from results/*.csv.   python3 scripts/plot.py --results results"""
import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import chi2

import plotstyle
from plotstyle import SERIES, MUTED, INK

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

ALGO_COLOR = {"dijkstra": SERIES[7], "astar": SERIES[0], "wastar_1.5": SERIES[3], "wastar_3": SERIES[4], "jps": SERIES[2]}
ALGO_LABEL = {"dijkstra": "Dijkstra", "astar": "A*", "wastar_1.5": "weighted A* (w=1.5)", "wastar_3": "weighted A* (w=3)", "jps": "JPS"}
TRACK_COLOR = {"pure_pursuit": SERIES[1], "stanley": SERIES[6], "lqr": SERIES[2]}


def pgm(path):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"P5"
        w, h = map(int, f.readline().split())
        f.readline()
        return np.frombuffer(f.read(), np.uint8).reshape(h, w)


def md(df, fmt=".2f"):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(f"{v:{fmt}}" if isinstance(v, (float, np.floating)) else str(v) for v in r) + " |")
    return "\n".join(lines)


def search(res):
    d = pd.read_csv(res / "search_compare.csv")
    order = ["dijkstra", "astar", "wastar_1.5", "wastar_3", "jps"]
    rows = []
    for size in sorted(d["size"].unique()):
        for a in order:
            x = d[(d["size"] == size) & (d.algo == a)]
            rows.append({"map": f"{size}x{size}", "algorithm": ALGO_LABEL[a], "median expansions": x.expansions.median(),
                         "median time [ms]": x.time_ms.median(), "mean cost / optimal": x.cost_ratio.mean(),
                         "worst cost / optimal": x.cost_ratio.max()})
    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.6))
    for i, (col, title) in enumerate([("expansions", "nodes expanded"), ("time_ms", "time [ms]")]):
        data = [d[(d["size"] == 512) & (d.algo == a)][col].to_numpy() for a in order]
        bp = ax[i].boxplot(data, orientation="horizontal", widths=0.55, patch_artist=True,
                           flierprops={"markersize": 3})
        for p, a in zip(bp["boxes"], order):
            p.set_facecolor(ALGO_COLOR[a])
            p.set_alpha(0.8)
        for m in bp["medians"]:
            m.set_color(INK)
        ax[i].set_xscale("log")
        ax[i].set_title(title)
        ax[i].set_yticks(range(1, 6), [ALGO_LABEL[a] for a in order] if i == 0 else [""] * 5)
    fig.suptitle("30 random 512x512 maps, random start and goal", fontsize=10)
    fig.tight_layout()
    fig.savefig(res / "search_compare.png")
    plt.close(fig)
    return md(pd.DataFrame(rows), ".3g")


def dstar(res):
    d = pd.read_csv(res / "dstar.csv")
    g = d.groupby("size")[["replans", "dstar_expansions", "astar_expansions", "dstar_ms", "astar_ms", "mismatch"]].sum()
    g["expansion ratio (A*/D*)"] = g.astar_expansions / g.dstar_expansions
    g["time ratio (A*/D*)"] = g.astar_ms / g.dstar_ms
    g = g.reset_index().rename(columns={"size": "map side", "dstar_expansions": "D* Lite expansions",
                                        "astar_expansions": "A* replan expansions", "dstar_ms": "D* Lite ms",
                                        "astar_ms": "A* ms", "mismatch": "cost mismatches"})
    truth = pgm(res / "dstar_truth.pgm")
    belief = pgm(res / "dstar_belief.pgm")
    tr = pd.read_csv(res / "dstar_trace.csv")
    ep = pd.read_csv(res / "dstar_endpoints.csv").iloc[0]
    h = truth.shape[0]
    img = np.full(truth.shape + (3,), 252, np.uint8)
    unknown_obst = (truth == 0) & (belief != 0)
    img[unknown_obst] = (200, 199, 192)  # never sensed
    img[belief == 0] = (31, 31, 30)      # discovered
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.imshow(img, extent=[0, truth.shape[1], 0, h], interpolation="nearest")
    ax.plot(tr.x + 0.5, tr.y + 0.5, color=SERIES[1], lw=1.6, label="robot path")
    ax.plot(ep.sx + 0.5, ep.sy + 0.5, "o", color=SERIES[0], ms=8, label="start")
    ax.plot(ep.gx + 0.5, ep.gy + 0.5, "*", color=SERIES[2], ms=13, label="goal")
    ax.set_title("D* Lite in an unknown map\n(dark: obstacles sensed, grey: never seen)", fontsize=10)
    ax.legend(loc="lower right", fontsize=8, frameon=True)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(res / "dstar_trace.png")
    plt.close(fig)
    return md(g, ".3g")


def hybrid(res):
    lot = pgm(res / "parking_lot.pgm")
    p = pd.read_csv(res / "hybrid_path.csv")
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    ax.imshow(np.where(lot == 0, 40, 252), cmap="gray", vmin=0, vmax=255, extent=[0, 60, 0, 40], interpolation="nearest")
    L, W, ro = 4.6, 1.9, 0.9
    for i in range(0, len(p), 25):
        r = p.iloc[i]
        c, s = np.cos(r.th), np.sin(r.th)
        corners = [(-ro, -W / 2), (L - ro, -W / 2), (L - ro, W / 2), (-ro, W / 2)]
        pts = [(r.x + c * a - s * b, r.y + s * a + c * b) for a, b in corners]
        ax.add_patch(Polygon(pts, closed=True, fill=False, lw=0.6,
                             edgecolor=SERIES[0] if r.dir > 0 else SERIES[1]))
    for sgn, col, lab in [(1, SERIES[0], "forward"), (-1, SERIES[1], "reverse")]:
        x = p[p.dir == sgn]
        ax.plot(x.x, x.y, ".", ms=1.5, color=col, label=lab)
    ax.set_title("Hybrid A*: rear-axle path and footprints, 0.1 m map", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, frameon=True, markerscale=6)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(res / "hybrid_parking.png")
    plt.close(fig)
    m = pd.read_csv(res / "hybrid_metrics.csv")
    return md(m, ".1f")


def rrt(res):
    h = pd.read_csv(res / "rrt_history.csv")
    ref = pd.read_csv(res / "rrt_reference.csv").astar_grid_cost_m.iloc[0]
    iters = np.arange(0, 6000, 25)
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.8), gridspec_kw={"width_ratios": [1.25, 1]})
    rows = []
    for v, col, lab in [("rrt_star", SERIES[1], "RRT*"), ("informed_rrt_star", SERIES[2], "Informed RRT*")]:
        curves = []
        for sd in range(1, 21):
            x = h[(h.variant == v) & (h.seed == sd)].sort_values("iter")
            c = np.full(len(iters), np.nan)
            for j, it in enumerate(iters):
                y = x[x.iter <= it]
                if len(y):
                    c[j] = y.cost.iloc[-1]
            curves.append(c)
        C = np.array(curves)
        med = np.nanmedian(C, axis=0)
        q1, q3 = np.nanpercentile(C, 25, axis=0), np.nanpercentile(C, 75, axis=0)
        ax[0].plot(iters, med, color=col, label=lab)
        ax[0].fill_between(iters, q1, q3, color=col, alpha=0.15, lw=0)
        for it in (1000, 2000, 4000, 5975):
            j = np.searchsorted(iters, it)
            rows.append({"planner": lab, "iterations": it, "median cost [m]": np.nanmedian(C[:, j]),
                         "solved [of 20]": int(np.sum(~np.isnan(C[:, j])))})
    ax[0].axhline(ref, color=MUTED, ls="--", lw=1)
    ax[0].annotate(f"8-connected A*: {ref:.1f} m", (3000, ref), xytext=(0, 4), textcoords="offset points",
                   fontsize=8, color=MUTED)
    ax[0].set_ylim(130, 160)
    ax[0].set_xlabel("iteration")
    ax[0].set_ylabel("best path cost [m]")
    ax[0].set_title("Median and IQR over 20 seeds")
    ax[0].legend(fontsize=8)
    cs = pgm(res / "rrt_cspace.pgm")
    ax[1].imshow(np.where(cs == 0, 40, 252), cmap="gray", vmin=0, vmax=255, extent=[0, 100, 0, 100], interpolation="nearest")
    for v, col, lab in [("rrt_star", SERIES[1], "RRT*"), ("informed_rrt_star", SERIES[2], "Informed RRT*")]:
        q = pd.read_csv(res / f"rrt_path_{v}.csv")
        ax[1].plot(q.x, q.y, color=col, lw=1.6, label=lab)
    ax[1].set_title("Final paths, seed 1 (inflated map)")
    ax[1].grid(False)
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    fig.tight_layout()
    fig.savefig(res / "rrt_convergence.png")
    plt.close(fig)
    return md(pd.DataFrame(rows), ".1f")


def tracking(res):
    d = pd.read_csv(res / "tracking_metrics.csv")
    g = d.groupby(["tracker", "state_source"]).agg(
        rms_lat_cm=("rms_lat_m", lambda v: 100 * v.mean()), max_lat_cm=("max_lat_m", lambda v: 100 * v.mean()),
        rms_head_deg=("rms_head_deg", "mean"), steer_rate=("steer_rate_rms", "mean"),
        est_rms_cm=("rms_pos_err_m", lambda v: 100 * v.mean())).reset_index()
    g.loc[g.state_source == "truth", "est_rms_cm"] = np.nan
    g.columns = ["tracker", "state used by controller", "RMS lateral [cm]", "max lateral [cm]",
                 "RMS heading [deg]", "RMS steer rate [rad/s]", "estimate RMS error [cm]"]

    road = pd.read_csv(res / "road.csv")
    fig, ax = plt.subplots(2, 1, figsize=(8.0, 6.4))
    ax[0].plot(road.x, road.y, color=MUTED, lw=4, alpha=0.35, label="reference")
    for tr in ["pure_pursuit", "stanley", "lqr"]:
        x = pd.read_csv(res / f"track_{tr}_ekf.csv")
        ax[0].plot(x.x, x.y, color=TRACK_COLOR[tr], lw=1.0, label=tr.replace("_", " "))
        ax[1].plot(x.t, 100 * x.e_lat, color=TRACK_COLOR[tr], lw=1.0, label=tr.replace("_", " "))
    x = pd.read_csv(res / "track_lqr_ekf_no_gate.csv")
    ax[1].plot(x.t, 100 * x.e_lat, color=SERIES[7], lw=0.9, alpha=0.8, label="lqr, GPS gate off")
    ax[0].set_aspect("equal")
    ax[0].set_title("640 m route with a lane change and a 40 m-radius curve, EKF state (seed 1)")
    ax[0].legend(fontsize=8, ncol=4, loc="upper left")
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[1].axvspan(24, 32, color="#e4e3dd", alpha=0.7, lw=0)
    ax[1].annotate("GPS dropout", (24.3, 85), fontsize=8, color=MUTED)
    ax[1].set_ylim(-100, 100)
    ax[1].set_xlabel("time [s]")
    ax[1].set_ylabel("true lateral error [cm]")
    ax[1].legend(fontsize=8, ncol=4, loc="lower left")
    fig.tight_layout()
    fig.savefig(res / "tracking.png")
    plt.close(fig)

    n = pd.read_csv(res / "ekf_nees.csv")
    N = int(n.runs.iloc[0])
    lo, hi = chi2.ppf(0.025, 3 * N) / N, chi2.ppf(0.975, 3 * N) / N
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    ax.plot(n.t, n.avg_nees, color=SERIES[0], lw=1.2, label=f"average NEES over {N} runs")
    ax.axhspan(lo, hi, color=SERIES[2], alpha=0.15, lw=0, label="95% band for a consistent filter")
    ax.axhline(3.0, color=MUTED, lw=0.8, ls=":")
    ax.axvspan(24, 32, color="#e4e3dd", alpha=0.7, lw=0)
    ax.annotate("GPS dropout", (24.3, ax.get_ylim()[1] * 0.92), fontsize=8, color=MUTED)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("NEES (x, y, heading)")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(res / "ekf_nees.png")
    plt.close(fig)
    inside = float(((n.avg_nees > lo) & (n.avg_nees < hi)).mean())
    nis = pd.read_csv(res / "ekf_nis.csv").nis
    phases = []
    for a, b, lab in [(2, 24, "before dropout"), (24, 32, "during dropout"), (32, 40, "8 s after reacquisition"), (40, 99, "rest")]:
        x = n[(n.t >= a) & (n.t < b)]
        phases.append({"phase": lab, "mean NEES": x.avg_nees.mean(),
                       "epochs inside 95% band": f"{100 * ((x.avg_nees > lo) & (x.avg_nees < hi)).mean():.0f}%"})
    summary = (f"Mean NIS over {len(nis)} accepted fixes: {nis.mean():.2f} (2.00 expected). "
               f"Average NEES inside the 95% band at {100 * inside:.0f}% of epochs (band {lo:.2f} to {hi:.2f}).")
    return md(g, ".2f"), md(pd.DataFrame(phases), ".2f"), summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    res = pathlib.Path(ap.parse_args().results)
    parts = ["# Results (generated by scripts/plot.py)"]
    parts.append("## Optimal grid search\n\n" + search(res))
    parts.append("## Replanning in unknown maps\n\n" + dstar(res))
    parts.append("## Hybrid A*\n\n" + hybrid(res))
    parts.append("## RRT* vs Informed RRT*\n\n" + rrt(res))
    t, ph, s = tracking(res)
    parts.append("## Tracking\n\n" + t + "\n\n" + s + "\n\n" + ph)
    (res / "summary.md").write_text("\n\n".join(parts) + "\n")
    print((res / "summary.md").read_text())


if __name__ == "__main__":
    main()
