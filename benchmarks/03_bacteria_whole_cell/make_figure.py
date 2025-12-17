#!/usr/bin/env python3
"""
Benchmark 03 - bacteria whole cell

Inputs:
  results/csvs/results_*.csv

Outputs:
  results/figures/
    - benchmark_n_cells.jpg
    - benchmark_cc3d_steps.jpg
    - benchmark_whole_cell_steps.jpg
    - benchmark_step_to_stationary_stage.jpg

Run: python make_figure.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Paths
SCRIPT_DIR = Path(__file__).parent
CSV_DIR = SCRIPT_DIR / "results" / "csvs"
OUT_DIR = SCRIPT_DIR / "results" / "figures"

NON_PARAM_COLS = {
    "execution_time",
    "execution_time_seconds",
    "time_reference",
    "outputdir",
    "trna_data_path",
    "mrna_data_path",
}


def _parse_execution_time_to_seconds(series):
    td = pd.to_timedelta(series, errors="coerce")
    return td.dt.total_seconds()


def _format_seconds(seconds):
    seconds_int = int(round(seconds))
    h = seconds_int // 3600
    m = (seconds_int % 3600) // 60
    s = seconds_int % 60
    if h > 0:
        return f"{h}h {m}m {s}s"
    return f"{m}m {s}s"


def _default_config(df):
    param_cols = [c for c in df.columns if c not in NON_PARAM_COLS]
    return df[param_cols].mode().iloc[0]


def _extract_isolated_curve(df, x_param):
    """
    Keep only rows where all other parameters match default config.
    """
    default = _default_config(df)
    other_params = [p for p in default.index.tolist() if p != x_param]
    mask = (df[other_params] == default[other_params]).all(axis=1)
    out = df[mask].copy()
    out = out.sort_values(by=x_param)
    return out


def _plot_param_vs_time(sources, x_param, out_path, title):
    fig, ax = plt.subplots(figsize=(12, 6))

    for source_name, df in sources.items():
        curve = _extract_isolated_curve(df, x_param)
        if curve.empty or curve[x_param].nunique() < 2:
            continue
        ax.plot(
            curve[x_param],
            curve["execution_time_seconds"],
            marker="o",
            linewidth=2,
            label=source_name,
        )
        for _, row in curve.iterrows():
            ax.text(
                row[x_param],
                row["execution_time_seconds"],
                f" {_format_seconds(row['execution_time_seconds'])}",
                fontsize=9,
                ha="left",
                va="center",
            )

    ax.set_title(title)
    ax.set_xlabel(x_param)
    ax.set_ylabel("Execution Time (seconds)")
    ax.grid(True, alpha=0.3)
    ax.legend(title="source")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_figures():
    csv_paths = sorted(CSV_DIR.glob("results_*.csv"))
    per_source = {}
    for p in csv_paths:
        df = pd.read_csv(p)
        df["execution_time_seconds"] = _parse_execution_time_to_seconds(df["execution_time"])
        df = df.dropna(subset=["execution_time_seconds"]).copy()
        per_source[p.stem] = df

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    outputs = []

    plots = [
        ("n_cells", OUT_DIR / "benchmark_n_cells.jpg", 'Impact of "n_cells" on Execution Time'),
        ("cc3d_steps", OUT_DIR / "benchmark_cc3d_steps.jpg", 'Impact of "cc3d_steps" on Execution Time'),
        ("whole_cell_steps", OUT_DIR / "benchmark_whole_cell_steps.jpg", 'Impact of "whole_cell_steps" on Execution Time'),
        ("step_to_stationary_stage", OUT_DIR / "benchmark_step_to_stationary_stage.jpg", 'Impact of "step_to_stationary_stage" on Execution Time'),
    ]

    for x_param, out_path, title in plots:
        _plot_param_vs_time(per_source, x_param, out_path, title)
        outputs.append(out_path)

    return outputs


def main():
    csv_paths = list(CSV_DIR.glob("results_*.csv"))
    if not csv_paths:
        print(f"No results_*.csv found in {CSV_DIR}")
        return 1

    outs = make_figures()
    for p in outs:
        print(f"Saved: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
