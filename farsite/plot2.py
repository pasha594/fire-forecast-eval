"""Pilot charts from runs/hindcast_scores.csv -> runs/pilot_skill.png
Row 1: IoU of predicted vs observed new growth at every perimeter, at its exact lead time
       (dots = individual run x perimeter, p50; line = 6-hour bin mean).
Row 2: predicted / observed new acres (pooled per 6-hour bin, log scale).
Row 3: satellite-pass capture (share of new detections inside the predicted area by then).
"""
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import RUNS

COL = {"pyrecast": "#2563eb", "farsite": "#dc2626"}
LAB = {"pyrecast": "PyreCast (ELMFIRE)", "farsite": "FARSITE (ours, automated)"}
BIN = 6


def paired(rows, truth, pct="50"):
    by = defaultdict(dict)
    for r in rows:
        if r["truth"] == truth and r["pct"] == pct and r.get("baseline", "common") == "common":
            by[(r["slug"], r["run_ts"], r["truth_id"], r["lead_h"])][r["model"]] = r
    return [v for v in by.values() if len(v) == 2 and all(x["passes"] == "1" for x in v.values())]


def main():
    rows = list(csv.DictReader(open(RUNS / "hindcast_scores.csv")))
    fires = sorted({r["slug"] for r in rows})
    fig, axes = plt.subplots(3, len(fires), figsize=(6 * len(fires), 10), constrained_layout=True, squeeze=False)
    for j, slug in enumerate(fires):
        sub = [r for r in rows if r["slug"] == slug]
        pp, hp = paired(sub, "perim"), paired(sub, "hotspot")
        edges = np.arange(0, 96 + BIN, BIN)
        mid = edges[:-1] + BIN / 2
        for m in ("pyrecast", "farsite"):
            lead = np.array([float(v[m]["lead_h"]) for v in pp])
            iou = np.array([float(v[m]["iou"] or 0) for v in pp])
            pred = np.array([float(v[m]["pred_new_ac"]) for v in pp])
            act = np.array([float(v[m]["act_new_ac"]) for v in pp])
            axes[0, j].scatter(lead, iou, s=8, alpha=0.25, color=COL[m])
            bi = np.digitize(lead, edges) - 1
            mean_iou = [iou[bi == k].mean() if (bi == k).any() else np.nan for k in range(len(mid))]
            over = [pred[bi == k].sum() / act[bi == k].sum() if (bi == k).any() and act[bi == k].sum() else np.nan
                    for k in range(len(mid))]
            axes[0, j].plot(mid, mean_iou, "-o", ms=3, color=COL[m], label=LAB[m])
            axes[1, j].plot(mid, over, "-o", ms=3, color=COL[m], label=LAB[m])
            hl = np.array([float(v[m]["lead_h"]) for v in hp])
            cap = np.array([float(v[m]["capture"]) for v in hp])
            hb = np.digitize(hl, edges) - 1
            axes[2, j].plot(mid, [cap[hb == k].mean() if (hb == k).any() else np.nan for k in range(len(mid))],
                            "-o", ms=3, color=COL[m], label=LAB[m])
        n_runs = len({v["pyrecast"]["run_ts"] for v in pp})
        axes[0, j].set_title(f"{slug}: IoU at each perimeter's exact lead (p50; {len(pp)} pairs, {n_runs} runs)", fontsize=10)
        axes[1, j].set_title(f"{slug}: predicted / observed new acres (p50, pooled per {BIN} h)", fontsize=10)
        axes[2, j].set_title(f"{slug}: new satellite detections inside forecast (p50; {len(hp)} passes)", fontsize=10)
        axes[1, j].axhline(1, color="#9ca3af", lw=0.8)
        axes[1, j].set_yscale("log")
        axes[0, j].set_ylim(0, 1)
        axes[2, j].set_ylim(0, 1)
        for ax in axes[:, j]:
            ax.set_xlim(0, 96)
            ax.set_xticks(range(0, 97, 12))
            ax.grid(alpha=0.25)
        axes[2, j].set_xlabel("lead time (hours after the run's start perimeter)")
    axes[0, 0].legend(fontsize=8, loc="upper right")
    fig.suptitle("FARSITE pilot vs PyreCast, same runs, scored at every observation time", fontsize=12)
    fig.savefig(RUNS / "pilot_skill.png", dpi=130)
    print(RUNS / "pilot_skill.png")


if __name__ == "__main__":
    main()
