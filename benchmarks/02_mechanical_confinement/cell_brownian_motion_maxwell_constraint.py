"""
Cell Brownian Motion with MaxwellMedium Constraint.

This script simulates cell motion under MaxwellMedium viscoelastic constraints:
- Volume and Contact plugins for shape control
- Directed Brownian motion: random direction changes every N MCS
- Records trajectory (trajectories.csv) for each run

Usage:
    python cell_brownian_motion_maxwell_constraint.py --single --elastic 0.5 --tau 1000
    python cell_brownian_motion_maxwell_constraint.py --sweep
    python cell_brownian_motion_maxwell_constraint.py --single --elastic 0.5 --tau 1000 --snapshots
"""

import os
import sys
import argparse
import numpy as np
import shutil
import subprocess
import multiprocessing as mp
import pandas as pd

from cc3d import CompuCellSetup
from cc3d.core.XMLUtils import ElementCC3D
from cc3d.core.PySteppables import SteppableBasePy


# ---------------------------
# Constants
# ---------------------------
FORCE_AMPLITUDE = 6.5  # Fixed amplitude for all experiments
NUM_RUNS = 100  # Number of runs per condition for statistics
NUM_WORKERS = 10  # Number of parallel workers
DIRECTION_CHANGE_INTERVAL = 100  # MCS between direction changes
SNAPSHOT_FREQUENCY = 100  # MCS between snapshots
LATTICE_SIZE = 96  # Potts lattice dimensions (cubic)


# ---------------------------
# XML-like configuration
# ---------------------------
def configure_simulation(
    number_of_processors=1,
    debug_output_frequency=10,
    potts_dimensions={"x": 96, "y": 96, "z": 96},
    potts_steps=10000,
    potts_temperature=10.0,
    neighbor_order=1,
    # Volume plugin parameters
    target_volume=2000.0,
    lambda_volume=50.0,
    # Contact energies (J=10 for all)
    contact_energy=10.0,
    # MaxwellMedium parameters
    maxwell_lambda_elastic=130.0,
    maxwell_reference_volume=10000.0,
    maxwell_relaxation_time=0.0,
):
    """Configure the CC3D simulation with Volume, Contact, and MaxwellMedium plugins."""
    cc3d = ElementCC3D("CompuCell3D", {"Version": "4.6.0"})

    # Metadata
    meta = cc3d.ElementCC3D("Metadata")
    meta.ElementCC3D("NumberOfProcessors", {}, number_of_processors)
    meta.ElementCC3D("DebugOutputFrequency", {}, debug_output_frequency)

    # Potts
    potts = cc3d.ElementCC3D("Potts")
    potts.ElementCC3D("Dimensions", potts_dimensions)
    potts.ElementCC3D("Steps", {}, potts_steps)
    potts.ElementCC3D("Temperature", {}, potts_temperature)
    potts.ElementCC3D("NeighborOrder", {}, neighbor_order)

    # CellType Plugin
    celltype = cc3d.ElementCC3D("Plugin", {"Name": "CellType"})
    celltype.ElementCC3D("CellType", {"TypeId": 0, "TypeName": "Medium"})
    celltype.ElementCC3D("CellType", {"TypeId": 1, "TypeName": "Cell"})

    # Volume Plugin
    volume = cc3d.ElementCC3D("Plugin", {"Name": "Volume"})
    volume.ElementCC3D("VolumeEnergyParameters", {
        "CellType": "Cell",
        "LambdaVolume": lambda_volume,
        "TargetVolume": target_volume
    })

    # Contact Plugin (J=10 for all interactions)
    contact = cc3d.ElementCC3D("Plugin", {"Name": "Contact"})
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Medium"}, contact_energy)
    contact.ElementCC3D("Energy", {"Type1": "Cell", "Type2": "Medium"}, contact_energy)
    contact.ElementCC3D("Energy", {"Type1": "Cell", "Type2": "Cell"}, contact_energy)
    contact.ElementCC3D("NeighborOrder", {}, neighbor_order)

    # ExternalPotential Plugin (for directed motion)
    cc3d.ElementCC3D("Plugin", {"Name": "ExternalPotential"})

    # CenterOfMass Plugin
    cc3d.ElementCC3D("Plugin", {"Name": "CenterOfMass"})

    # MaxwellMedium Plugin
    maxmed = cc3d.ElementCC3D("Plugin", {"Name": "MaxwellMedium"})
    maxmed.ElementCC3D("LambdaElastic", {}, maxwell_lambda_elastic)
    maxmed.ElementCC3D("ReferenceVolume", {}, maxwell_reference_volume)
    maxmed.ElementCC3D("TauR", {}, maxwell_relaxation_time)

    CompuCellSetup.setSimulationXMLDescription(cc3d)


# ---------------------------
# Steppables
# ---------------------------

class InitialLayoutSteppable(SteppableBasePy):
    """
    Creates a single spherical cell at the center of the lattice.
    """

    def __init__(self, frequency=1, radius=3):
        super().__init__(frequency)
        self.radius = radius

    def start(self):
        # Place cell at the exact center of the lattice
        cx = self.dim.x // 2
        cy = self.dim.y // 2
        cz = self.dim.z // 2

        cell = self.new_cell(self.CELL)

        # Create a spherical cell
        for dx in range(-self.radius, self.radius + 1):
            for dy in range(-self.radius, self.radius + 1):
                for dz in range(-self.radius, self.radius + 1):
                    if dx*dx + dy*dy + dz*dz <= self.radius**2:
                        x = max(0, min(self.dim.x - 1, cx + dx))
                        y = max(0, min(self.dim.y - 1, cy + dy))
                        z = max(0, min(self.dim.z - 1, cz + dz))
                        self.cell_field[x, y, z] = cell

        print(f"Cell created at center ({cx}, {cy}, {cz}) with radius {self.radius}")


class BrownianMotionSteppable(SteppableBasePy):
    """
    Implements directed Brownian motion:
    - Random direction in 3D, changed every N MCS
    - Force applied via ExternalPotential (lambdaVec)
    - Records trajectory to CSV
    """

    def __init__(
        self,
        frequency=1,
        force_amplitude=FORCE_AMPLITUDE,
        direction_change_interval=DIRECTION_CHANGE_INTERVAL,
        outputdir="./output/",
        record_frequency=10,
        snapshot_frequency=SNAPSHOT_FREQUENCY,
        save_snapshots=False,
    ):
        super().__init__(frequency)

        self.force_amplitude = force_amplitude
        self.direction_change_interval = direction_change_interval
        self.outputdir = outputdir
        self.record_frequency = record_frequency
        self.snapshot_frequency = snapshot_frequency
        self.save_snapshots = save_snapshots

        self.trajectories = []
        self.current_direction = None

        os.makedirs(self.outputdir, exist_ok=True)

    def _random_direction(self):
        """Generate a random unit vector in 3D (uniform on sphere)."""
        # Use Marsaglia's method for uniform distribution on sphere
        while True:
            x1 = np.random.uniform(-1, 1)
            x2 = np.random.uniform(-1, 1)
            s = x1**2 + x2**2
            if s < 1:
                break

        sqrt_term = np.sqrt(1 - s)
        dx = 2 * x1 * sqrt_term
        dy = 2 * x2 * sqrt_term
        dz = 1 - 2 * s

        return np.array([dx, dy, dz])

    def _save_instance_segmentation(self, mcs):
        """Save the current cell field state as a TIFF image."""
        from skimage.io import imsave

        Z, Y, X = self.dim.z, self.dim.y, self.dim.x
        img = np.zeros((Z, Y, X), dtype=np.uint16)

        cf = self.cell_field
        for z in range(Z):
            for y in range(Y):
                for x in range(X):
                    cell = cf[x, y, z]
                    if cell is not None:
                        img[z, y, x] = cell.id

        name = "snapshot_" + str(mcs).zfill(7) + ".tiff"
        imsave(os.path.join(self.outputdir, name), img.astype(np.uint8), check_contrast=False)

    def start(self):
        # Initialize with a random direction
        self.current_direction = self._random_direction()

        for cell in self.cell_list:
            # Apply initial force direction
            cell.lambdaVecX = -self.force_amplitude * self.current_direction[0]
            cell.lambdaVecY = -self.force_amplitude * self.current_direction[1]
            cell.lambdaVecZ = -self.force_amplitude * self.current_direction[2]

            # Record initial position
            self.trajectories.append({
                'mcs': 0,
                'cell_id': cell.id,
                'x': cell.xCOM,
                'y': cell.yCOM,
                'z': cell.zCOM,
                'dir_x': self.current_direction[0],
                'dir_y': self.current_direction[1],
                'dir_z': self.current_direction[2]
            })

            print(f"Initial direction: ({self.current_direction[0]:.3f}, "
                  f"{self.current_direction[1]:.3f}, {self.current_direction[2]:.3f})")

    def step(self, mcs: int):
        # Change direction every N MCS
        if mcs > 0 and mcs % self.direction_change_interval == 0:
            self.current_direction = self._random_direction()
            print(f"MCS {mcs}: New direction ({self.current_direction[0]:.3f}, "
                  f"{self.current_direction[1]:.3f}, {self.current_direction[2]:.3f})")

        for cell in self.cell_list:
            # Apply force in current direction (negative for attraction)
            cell.lambdaVecX = -self.force_amplitude * self.current_direction[0]
            cell.lambdaVecY = -self.force_amplitude * self.current_direction[1]
            cell.lambdaVecZ = -self.force_amplitude * self.current_direction[2]

            # Record trajectory
            if mcs % self.record_frequency == 0:
                self.trajectories.append({
                    'mcs': mcs,
                    'cell_id': cell.id,
                    'x': cell.xCOM,
                    'y': cell.yCOM,
                    'z': cell.zCOM,
                    'dir_x': self.current_direction[0],
                    'dir_y': self.current_direction[1],
                    'dir_z': self.current_direction[2]
                })

                if mcs % (self.record_frequency * 50) == 0:
                    print(f"MCS {mcs}: COM=({cell.xCOM:.2f}, {cell.yCOM:.2f}, {cell.zCOM:.2f})")

        # Save snapshot at regular intervals (only if enabled)
        if self.save_snapshots and mcs % self.snapshot_frequency == 0:
            self._save_instance_segmentation(mcs)

    def finish(self):
        df = pd.DataFrame(self.trajectories)
        trajectory_file = os.path.join(self.outputdir, "trajectories.csv")
        df.to_csv(trajectory_file, index=False)
        print(f"Trajectories saved to {trajectory_file}")


# ---------------------------
# Summary generation
# ---------------------------

def generate_trajectories_summary(data_dir: str, csv_output_path: str) -> str:
    """
    Consolidate all trajectories.csv files into a single trajectories_summary.csv.
    
    The summary contains columns: tau, elastic, run, mcs, x, y, z
    """
    import re

    tau_re = re.compile(r"tau_(\d+)")
    elastic_re = re.compile(r"elastic_([0-9.]+)")
    run_re = re.compile(r"run_(\d+)")

    rows = []
    data_path = os.path.abspath(data_dir)

    for root, dirs, files in os.walk(data_path):
        if "trajectories.csv" not in files:
            continue
        
        rel = os.path.relpath(root, data_path)
        
        tau_match = tau_re.search(rel)
        elastic_match = elastic_re.search(rel)
        run_match = run_re.search(rel)
        
        if not (tau_match and elastic_match and run_match):
            continue
        
        tau = int(tau_match.group(1))
        elastic = float(elastic_match.group(1))
        run_idx = int(run_match.group(1))
        
        traj_df = pd.read_csv(os.path.join(root, "trajectories.csv"))
        for _, row in traj_df.iterrows():
            rows.append({
                "tau": tau,
                "elastic": elastic,
                "run": run_idx,
                "mcs": int(row["mcs"]),
                "x": float(row["x"]),
                "y": float(row["y"]),
                "z": float(row["z"]),
            })

    summary_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(csv_output_path), exist_ok=True)
    summary_df.to_csv(csv_output_path, index=False)
    print(f"Trajectories summary saved to {csv_output_path} ({len(summary_df)} rows)")
    return csv_output_path


# ---------------------------
# Simulation Runners
# ---------------------------

def run_single_cc3d_simulation(
    elastic_modulus=0.0,
    relaxation_time=0.0,
    output_dir="./output",
    run_index=0,
    save_snapshots=False,
):
    """Run a single CC3D simulation (called via subprocess)."""
    run_dir = os.path.join(output_dir, f"run_{run_index}")
    os.makedirs(run_dir, exist_ok=True)

    potts_steps = 5000

    configure_simulation(
        number_of_processors=1,
        debug_output_frequency=100,
        potts_dimensions={"x": 96, "y": 96, "z": 96},
        potts_steps=potts_steps,
        potts_temperature=10.0,
        neighbor_order=1,
        target_volume=2000.0,
        lambda_volume=2.0,
        contact_energy=10.0,
        maxwell_lambda_elastic=elastic_modulus,
        maxwell_reference_volume=10000.0,
        maxwell_relaxation_time=relaxation_time,
    )

    # Register steppables
    CompuCellSetup.register_steppable(
        steppable=InitialLayoutSteppable(frequency=1, radius=3)
    )

    CompuCellSetup.register_steppable(
        steppable=BrownianMotionSteppable(
            frequency=1,
            force_amplitude=FORCE_AMPLITUDE,
            direction_change_interval=DIRECTION_CHANGE_INTERVAL,
            outputdir=run_dir,
            record_frequency=10,
            save_snapshots=save_snapshots,
        )
    )

    # Run simulation
    CompuCellSetup.run()

    return run_dir


def _run_single_worker(args):
    """
    Worker function for parallel execution of single runs.
    Args is a tuple: (run_idx, elastic_modulus, relaxation_time, output_dir, script_path, python_exe, save_snapshots)
    """
    run_idx, elastic_modulus, relaxation_time, output_dir, script_path, python_exe, save_snapshots = args
    
    print(f"  Starting run {run_idx} (elastic={elastic_modulus}, tau={relaxation_time})")
    
    # Launch simulation in subprocess
    cmd = [
        python_exe, script_path, '--run_single',
        '--elastic', str(elastic_modulus),
        '--tau', str(relaxation_time),
        '--output', output_dir,
        '--run_index', str(run_idx),
    ]
    if save_snapshots:
        cmd.append('--snapshots')
    
    result = subprocess.run(cmd, capture_output=True)
    
    if result.returncode != 0:
        print(f"  WARNING: Run {run_idx} failed!")
        return None
    
    run_dir = os.path.join(output_dir, f"run_{run_idx}")
    traj_file = os.path.join(run_dir, "trajectories.csv")
    
    if os.path.exists(traj_file):
        print(f"  Completed run {run_idx}")
        return traj_file
    
    return None


def run_multi_simulation(
    elastic_modulus=0.0,
    relaxation_time=0.0,
    output_base_dir="./single_test",
    num_runs=NUM_RUNS,
    num_workers=NUM_WORKERS,
    save_snapshots=False,
):
    """
    Run multiple simulations.
    Each run is executed in a separate subprocess, parallelized across workers.
    """
    # Create output directory
    output_dir = os.path.join(
        output_base_dir,
        f"tau_{relaxation_time:.0f}",
        f"elastic_{elastic_modulus:.1f}"
    )
    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    script_path = os.path.abspath(__file__)
    python_exe = sys.executable

    print(f"\n{'='*60}")
    print(f"Running {num_runs} simulations in parallel ({num_workers} workers)")
    print(f"elastic={elastic_modulus}, tau={relaxation_time}")
    print(f"{'='*60}\n")

    # Prepare arguments for parallel workers
    worker_args = [
        (run_idx, elastic_modulus, relaxation_time, output_dir, script_path, python_exe, save_snapshots)
        for run_idx in range(num_runs)
    ]

    # Run in parallel using multiprocessing Pool
    with mp.Pool(processes=num_workers) as pool:
        results = pool.map(_run_single_worker, worker_args)

    # Collect successful trajectory files
    trajectory_files = [f for f in results if f is not None]
    print(f"\nCompleted {len(trajectory_files)}/{num_runs} runs successfully")

    return output_dir


def run_parameter_sweep(output_base_dir="./results/sweep/data", save_snapshots=False):
    """
    Run a parameter sweep over elastic_modulus and relaxation_time.
    Each condition runs NUM_RUNS simulations.
    At the end, generates a consolidated trajectories_summary.csv.
    """
    if os.path.exists(output_base_dir):
        shutil.rmtree(output_base_dir)
    os.makedirs(output_base_dir, exist_ok=True)

    # Parameters to sweep
    elastic_moduli = [0, 0.33, 0.66, 1.0, 1.33]
    relaxation_times = [0, 1000, 5000, 10000]

    total_conditions = len(elastic_moduli) * len(relaxation_times)
    condition_count = 0

    for tau in relaxation_times:
        for elastic in elastic_moduli:
            condition_count += 1
            print(f"\n{'#'*60}")
            print(f"Condition {condition_count}/{total_conditions}")
            print(f"tau={tau}, elastic={elastic}")
            print(f"{'#'*60}\n")

            run_multi_simulation(
                elastic_modulus=elastic,
                relaxation_time=tau,
                output_base_dir=output_base_dir,
                num_runs=NUM_RUNS,
                save_snapshots=save_snapshots,
            )

    # Generate consolidated summary CSV
    csv_dir = os.path.join(os.path.dirname(output_base_dir), "csvs")
    summary_path = os.path.join(csv_dir, "trajectories_summary.csv")
    generate_trajectories_summary(output_base_dir, summary_path)

    print("\nParameter sweep completed!")


# ---------------------------
# Main
# ---------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Test MaxwellMedium plugin with Brownian directed motion"
    )
    parser.add_argument(
        '--sweep',
        action='store_true',
        help='Run parameter sweep (multiple conditions, each with NUM_RUNS simulations)'
    )
    parser.add_argument(
        '--single',
        action='store_true',
        help='Run multi-simulation for a single condition (NUM_RUNS runs)'
    )
    parser.add_argument(
        '--run_single',
        action='store_true',
        help='Internal: run a single CC3D simulation (used by subprocess)'
    )
    parser.add_argument(
        '--elastic',
        type=float,
        default=0.0,
        help='Elastic modulus (LambdaElastic) for MaxwellMedium (default: 0.0)'
    )
    parser.add_argument(
        '--tau',
        type=float,
        default=0.0,
        help='Relaxation time (TauR) for MaxwellMedium (default: 0.0)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='./results/sweep/data',
        help='Output base directory'
    )
    parser.add_argument(
        '--run_index',
        type=int,
        default=0,
        help='Internal: run index for subprocess'
    )
    parser.add_argument(
        '--num_runs',
        type=int,
        default=NUM_RUNS,
        help=f'Number of runs per condition (default: {NUM_RUNS})'
    )
    parser.add_argument(
        '--snapshots',
        action='store_true',
        help='Enable saving TIFF snapshots (disabled by default)'
    )

    args = parser.parse_args()

    if args.sweep:
        print("Running parameter sweep...")
        print(f"Each condition will run {NUM_RUNS} simulations")
        print(f"Snapshots: {'enabled' if args.snapshots else 'disabled'}")
        run_parameter_sweep(output_base_dir=args.output, save_snapshots=args.snapshots)
    elif args.run_single:
        # Internal: run single CC3D simulation (called by subprocess)
        run_single_cc3d_simulation(
            elastic_modulus=args.elastic,
            relaxation_time=args.tau,
            output_dir=args.output,
            run_index=args.run_index,
            save_snapshots=args.snapshots,
        )
    elif args.single:
        print(f"Running multi-simulation ({args.num_runs} runs):")
        print(f"  elastic_modulus = {args.elastic}")
        print(f"  tau = {args.tau}")
        print(f"  force_amplitude = {FORCE_AMPLITUDE} (fixed)")
        print(f"  snapshots: {'enabled' if args.snapshots else 'disabled'}")
        
        output_dir = run_multi_simulation(
            elastic_modulus=args.elastic,
            relaxation_time=args.tau,
            output_base_dir=args.output,
            num_runs=args.num_runs,
            save_snapshots=args.snapshots,
        )
        
        print(f"\nSimulation completed! Output in {output_dir}")
    else:
        # Default: run multi-simulation test
        print(f"Running multi-simulation test ({NUM_RUNS} runs):")
        print(f"  elastic_modulus = {args.elastic}")
        print(f"  tau = {args.tau}")
        print(f"  force_amplitude = {FORCE_AMPLITUDE} (fixed)")
        print(f"  snapshots: {'enabled' if args.snapshots else 'disabled'}")
        
        output_dir = run_multi_simulation(
            elastic_modulus=args.elastic,
            relaxation_time=args.tau,
            output_base_dir=args.output,
            num_runs=NUM_RUNS,
            save_snapshots=args.snapshots,
        )
        
        print(f"\nTest completed! Output in {output_dir}")
