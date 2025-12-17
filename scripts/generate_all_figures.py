#!/usr/bin/env python3
"""
Generate all benchmark figures.

Calls make_figure.py for each benchmark:
  - 01_spheroid_scalability
  - 02_mechanical_confinement
  - 03_bacteria_whole_cell

Run from repo root:
  python scripts/generate_all_figures.py
"""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
BENCHMARKS = [
    "01_spheroid_scalability",
    "02_mechanical_confinement",
    "03_bacteria_whole_cell",
]


def main():
    failed = []

    for benchmark in BENCHMARKS:
        script = REPO_ROOT / "benchmarks" / benchmark / "make_figure.py"
        if not script.exists():
            print(f"[SKIP] {benchmark}: make_figure.py not found")
            continue

        print(f"\n{'='*60}")
        print(f"[RUN] {benchmark}")
        print(f"{'='*60}")

        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=script.parent,
        )

        if result.returncode != 0:
            print(f"[FAIL] {benchmark}")
            failed.append(benchmark)
        else:
            print(f"[OK] {benchmark}")

    print(f"\n{'='*60}")
    if failed:
        print(f"Failed: {', '.join(failed)}")
        return 1
    print("All figures generated successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

