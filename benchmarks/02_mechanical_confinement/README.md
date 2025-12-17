# Benchmark 02 — Maxwell Medium Confinement

## Objective

This benchmark studies the **constrained motion** of a single cell under a **Maxwell viscoelastic medium**. The Maxwell model combines an elastic response with stress relaxation, providing a minimal proxy for **tumor cell confinement** in surrounding tissue or extracellular matrix.

## Method

A single spherical cell undergoes directed Brownian motion (random direction changes every 100 MCS) while subjected to the MaxwellMedium plugin. The simulation sweeps over two parameters:

- **τ (relaxation time):** 0, 1000, 5000, 10000 MCS
- **λ_elastic (elastic modulus):** 0, 0.33, 0.66, 1.0, 1.33

For each (τ, λ_elastic) condition, 100 independent runs are performed to gather statistics on cell trajectories.

## How to run

```bash
# Full parameter sweep
python cell_brownian_motion_maxwell_constraint.py --sweep

# Single condition (example)
python cell_brownian_motion_maxwell_constraint.py --single --elastic 0.5 --tau 1000
```

## Results

### 3D Trajectories

Each figure shows 100 superimposed trajectories for a given (τ, elastic) condition. Higher elasticity confines motion toward the center; higher relaxation time allows partial recovery of mobility.

| τ = 0 | τ = 1000 | τ = 5000 | τ = 10000 |
|:-----:|:--------:|:--------:|:---------:|
| ![](results/sweep/figures/trajectories_3d_tau_0.0_elastic_0.0.png) | ![](results/sweep/figures/trajectories_3d_tau_1000.0_elastic_0.0.png) | ![](results/sweep/figures/trajectories_3d_tau_5000.0_elastic_0.0.png) | ![](results/sweep/figures/trajectories_3d_tau_10000.0_elastic_0.0.png) |
| ![](results/sweep/figures/trajectories_3d_tau_0.0_elastic_0.7.png) | ![](results/sweep/figures/trajectories_3d_tau_1000.0_elastic_0.7.png) | ![](results/sweep/figures/trajectories_3d_tau_5000.0_elastic_0.7.png) | ![](results/sweep/figures/trajectories_3d_tau_10000.0_elastic_0.7.png) |
| ![](results/sweep/figures/trajectories_3d_tau_0.0_elastic_1.3.png) | ![](results/sweep/figures/trajectories_3d_tau_1000.0_elastic_1.3.png) | ![](results/sweep/figures/trajectories_3d_tau_5000.0_elastic_1.3.png) | ![](results/sweep/figures/trajectories_3d_tau_10000.0_elastic_1.3.png) |

*Rows: elastic = 0.0, 0.7, 1.3 (top to bottom). Columns: τ = 0, 1000, 5000, 10000.*

### Radial distribution

The radial distribution p(r) measures how far from the lattice center cell positions are observed. Confinement (high elasticity, low τ) shifts the distribution toward smaller radii.

| τ = 0 | τ = 1000 |
|:-----:|:--------:|
| ![](results/sweep/figures/radial_distribution_tau_0.png) | ![](results/sweep/figures/radial_distribution_tau_1000.png) |

| τ = 5000 | τ = 10000 |
|:--------:|:---------:|
| ![](results/sweep/figures/radial_distribution_tau_5000.png) | ![](results/sweep/figures/radial_distribution_tau_10000.png) |

## Output files

- `results/sweep/csvs/trajectories_summary.csv` — consolidated trajectory data
- `results/sweep/figures/` — all figures above

