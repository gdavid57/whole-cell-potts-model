#!/usr/bin/env python3
"""
Benchmark 02 - mechanical confinement

Generates figures from trajectories_summary.csv:
  1. Radial distribution of coordinates (one figure per tau)
  2. 3D trajectory plots (one figure per tau/elastic experiment)

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
SUMMARY_CSV = SCRIPT_DIR / "results" / "sweep" / "csvs" / "trajectories_summary.csv"
OUT_DIR = SCRIPT_DIR / "results" / "sweep" / "figures"
LATTICE_SIZE = 96
N_BINS = 60


def make_radial_distributions(df):
    """
    Generate radial distribution figures.
    One figure per tau value, with one curve per elastic value.
    """
    center = LATTICE_SIZE / 2.0
    r_max = np.sqrt(3.0) * center
    bins = np.linspace(0.0, r_max, N_BINS + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_widths = np.diff(bins)

    df = df.copy()
    df["r"] = np.sqrt(
        (df["x"] - center) ** 2 +
        (df["y"] - center) ** 2 +
        (df["z"] - center) ** 2
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    saved = []

    for tau in sorted(df["tau"].unique()):
        fig, ax = plt.subplots(figsize=(10, 6))
        df_tau = df[df["tau"] == tau]

        for elastic in sorted(df_tau["elastic"].unique()):
            df_cond = df_tau[df_tau["elastic"] == elastic]
            counts, _ = np.histogram(df_cond["r"].values, bins=bins)
            total = counts.sum()
            if total == 0:
                continue
            pdf = (counts / total) / bin_widths
            ax.plot(bin_centers, pdf, linewidth=2, label=f"elastic={elastic:g}")

        ax.set_xlabel("Radius r (pixels)")
        ax.set_ylabel("Radial distribution p(r) (1/pixels)")
        ax.set_title(f"Radial distribution of positions (tau={tau})")
        ax.grid(True, alpha=0.3)
        ax.legend(title="Condition", fontsize=9)

        out_path = OUT_DIR / f"radial_distribution_tau_{tau}.png"
        fig.tight_layout()
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_path)

    return saved


def make_trajectories_3d(df):
    """
    Generate 3D trajectory plots for each (tau, elastic) experiment.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    saved = []

    conditions = df.groupby(["tau", "elastic"]).size().reset_index()[["tau", "elastic"]]

    for _, row in conditions.iterrows():
        tau = row["tau"]
        elastic = row["elastic"]

        df_cond = df[(df["tau"] == tau) & (df["elastic"] == elastic)]
        if df_cond.empty:
            continue

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")

        runs = df_cond["run"].unique()
        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]

        for i, run_id in enumerate(runs):
            df_run = df_cond[df_cond["run"] == run_id].sort_values("mcs")
            ax.plot(
                df_run["x"].values,
                df_run["y"].values,
                df_run["z"].values,
                linewidth=0.5,
                alpha=0.6,
                color=colors[i % len(colors)],
            )

        ax.set_xlim(0, LATTICE_SIZE)
        ax.set_ylim(0, LATTICE_SIZE)
        ax.set_zlim(0, LATTICE_SIZE)

        ax.set_xlabel("X (pixels)")
        ax.set_ylabel("Y (pixels)")
        ax.set_zlabel("Z (pixels)")
        ax.set_title(f"3D Trajectories (tau={tau}, elastic={elastic:g})\n{len(runs)} runs")

        out_path = OUT_DIR / f"trajectories_3d_tau_{tau}_elastic_{elastic}.png"
        fig.tight_layout()
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_path)

    return saved


def main():
    if not SUMMARY_CSV.exists():
        print(f"Summary CSV not found: {SUMMARY_CSV}")
        return 1

    print(f"Loading {SUMMARY_CSV}...")
    df = pd.read_csv(SUMMARY_CSV)
    if df.empty:
        print("Summary CSV is empty.")
        return 1
    print(f"  {len(df)} rows loaded.")

    print("\nGenerating radial distribution figures...")
    for p in make_radial_distributions(df):
        print(f"  Saved: {p}")

    print("\nGenerating 3D trajectory figures...")
    for p in make_trajectories_3d(df):
        print(f"  Saved: {p}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
