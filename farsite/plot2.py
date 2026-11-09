"""Pilot charts from runs/hindcast_scores.csv (score2.py) -> runs/pilot_skill.png
Rows matched across every model shown (same run, percentile 50, perimeter), each at its exact
lead (hours after the run): dots are individual scores, lines join each 24 h day's mean.
Row 1: IoU of predicted vs observed new growth.
Row 2: predicted / observed new acres, pooled per day (log scale).
Row 3: hotspot capture (share of the window's new detections inside the forecast).
"""
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import RUNS

COL = {"pyrecast": "#2563eb", "farsite": "#dc2626", "farsite-cal": "#16a34a", "farsite-calr": "#9333ea",
       "farsite-calbar": "#ca8a04", "farsite-calbarr": "#0891b2"}
LAB = {"pyrecast": "PyreCast (ELMFIRE)", "farsite": "FARSITE v1 (gusts, uncalibrated)",
       "farsite-cal": "FARSITE calibrated (constant factor)",
       "farsite-calr": "FARSITE calibrated (factor eases off)",
       "farsite-calbar": "FARSITE calibrated + control lines (constant)",
       "farsite-calbarr": "FARSITE calibrated + control lines (eases off)"}
f = lambda x: float(x) if x not in ("", None) else np.nan


def matched(rows, models, pct="50"):
    by = defaultdict(dict)
    for r in rows:
        if r["percentile"] == pct:
            by[(r["slug"], r["run_ts"], r["actual_perim_ms"])][r["model"]] = r
    return [v for v in by.values() if all(m in v for m in models)]


def main():
    rows = list(csv.DictReader(open(RUNS / "hindcast_scores.csv")))
    fires = sorted({r["slug"] for r in rows})
    fig, axes = plt.subplots(3, len(fires), figsize=(6.5 * len(fires), 10), constrained_layout=True, squeeze=False)
    for j, slug in enumerate(fires):
        sub = [r for r in rows if r["slug"] == slug]
        models = [m for m in COL if any(r["model"] == m for r in sub)]
        mm = matched(sub, models)
        days = sorted({int(v[models[0]]["day"]) for v in mm})
        for m in models:
            lead = np.array([f(v[m]["lead_h"]) for v in mm])
            iou = np.array([f(v[m]["iou"]) for v in mm])
            cap = np.array([f(v[m]["hs_capture"]) for v in mm])
            pred = np.array([f(v[m]["pred_new_acres"]) for v in mm])
            act = np.array([f(v[m]["act_new_acres"]) for v in mm])
            day = np.array([int(v[m]["day"]) for v in mm])
            mid = [(d - 0.5) * 24 for d in days]
            axes[0, j].scatter(lead, iou, s=7, alpha=0.25, color=COL[m])
            axes[0, j].plot(mid, [np.nanmean(iou[day == d]) for d in days], "-o", ms=3, color=COL[m], label=LAB[m])
            axes[1, j].plot(mid, [pred[day == d].sum() / act[day == d].sum() if act[day == d].sum() else np.nan for d in days],
                            "-o", ms=3, color=COL[m], label=LAB[m])
            axes[2, j].scatter(lead, cap, s=7, alpha=0.2, color=COL[m])
            axes[2, j].plot(mid, [np.nanmean(cap[day == d]) for d in days], "-o", ms=3, color=COL[m], label=LAB[m])
        n_runs = len({v[models[0]]["run_ts"] for v in mm})
        axes[0, j].set_title(f"{slug}: IoU at each perimeter's own time (p50; {len(mm)} perimeters, {n_runs} runs)", fontsize=10)
        axes[1, j].set_title(f"{slug}: predicted / observed new acres (p50, pooled per day)", fontsize=10)
        axes[2, j].set_title(f"{slug}: new satellite detections inside the forecast (p50)", fontsize=10)
        axes[1, j].axhline(1, color="#9ca3af", lw=0.8)
        axes[1, j].set_yscale("log")
        top = np.nanmax([f(v[m]["iou"]) for v in mm for m in models] or [1])
        axes[0, j].set_ylim(0, max(0.05, min(1, 1.15 * top)))
        axes[2, j].set_ylim(0, 1)
        xmax = max(168, max((f(v[models[0]]["lead_h"]) for v in mm), default=168))
        for ax in axes[:, j]:
            ax.set_xlim(0, xmax)
            ax.set_xticks(range(0, int(xmax) + 1, 24))
            ax.set_xticklabels([f"{h // 24}d" if h else "run" for h in range(0, int(xmax) + 1, 24)])
            ax.grid(alpha=0.25)
        axes[2, j].set_xlabel("time after the run (each dot at its perimeter's exact time)")
    axes[0, 0].legend(fontsize=8, loc="upper right")
    fig.suptitle("FARSITE pilot vs PyreCast, scored the dashboard's way (same perimeters, one common start)", fontsize=12)
    fig.savefig(RUNS / "pilot_skill.png", dpi=130)
    print(RUNS / "pilot_skill.png")


if __name__ == "__main__":
    main()
