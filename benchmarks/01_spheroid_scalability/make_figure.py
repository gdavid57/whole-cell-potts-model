#!/usr/bin/env python3
"""
Benchmark 01 - spheroid scalability

Generates figures from:
  results/csvs/generation_times.csv

Outputs to:
  results/figures/scalability_generation_times.png

Run: python make_figure.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Paths
SCRIPT_DIR = Path(__file__).parent
CSV_PATH = SCRIPT_DIR / "results" / "csvs" / "generation_times.csv"
OUT_DIR = SCRIPT_DIR / "results" / "figures"

# Hardware
NUM_LOGICAL_CORES = 128  # AMD CPU used for this benchmark


def make_figures():
    df = pd.read_csv(CSV_PATH)
    required = {"cell_count", "elapsed_seconds"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {CSV_PATH}: {sorted(missing)}")

    df = df.sort_values("cell_count").reset_index(drop=True)
    df["delta_elapsed_seconds"] = df["elapsed_seconds"].diff()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "scalability_generation_times.png"

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Execution time vs number of cells
    ax = axes[0]
    ax.plot(df["cell_count"], df["elapsed_seconds"], marker="o", linewidth=2)
    ax.axvline(NUM_LOGICAL_CORES, color="red", linestyle="--", linewidth=0.8, label=f"{NUM_LOGICAL_CORES} cores")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Execution time (seconds)")
    ax.set_title("Execution time vs number of cells")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Delta execution time vs number of cells
    ax = axes[1]
    ax.plot(df["cell_count"], df["delta_elapsed_seconds"], marker="o", linewidth=2)
    ax.axvline(NUM_LOGICAL_CORES, color="red", linestyle="--", linewidth=0.8, label=f"{NUM_LOGICAL_CORES} cores")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Number of cells")
    ax.set_ylabel(r"$\Delta$ execution time (seconds)")
    ax.set_title(r"$\Delta$ execution time vs number of cells")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # annotate points
    for _, row in df.iterrows():
        ax0_y = row["elapsed_seconds"]
        axes[0].annotate(
            f"{ax0_y:.0f}s",
            (row["cell_count"], ax0_y),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
        )
        if np.isfinite(row["delta_elapsed_seconds"]):
            ax1_y = row["delta_elapsed_seconds"]
            axes[1].annotate(
                f"{ax1_y:.0f}s",
                (row["cell_count"], ax1_y),
                textcoords="offset points",
                xytext=(6, 6),
                fontsize=8,
            )

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return out_path


def main():
    if not CSV_PATH.exists():
        print(f"CSV not found: {CSV_PATH}")
        return 1

    out = make_figures()
    print(f"Saved: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
