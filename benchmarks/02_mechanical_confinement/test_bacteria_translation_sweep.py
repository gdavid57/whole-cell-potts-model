"""
MaxwellMedium Plugin Test with Brownian Directed Motion and MSD Analysis.

This script tests the MaxwellMedium plugin using a single cell with:
- Volume and Contact plugins for shape control (no BacteriaShape)
- Directed Brownian motion: random direction changes every N MCS
- Mean Square Displacement (MSD) and effective diffusion coefficient measurement
- Multiple runs (M=5) per condition for statistical analysis

Usage:
    python test_bacteria_translation_sweep.py --single --elastic 100.0 --tau 1000.0
    python test_bacteria_translation_sweep.py --sweep
"""

import os
import sys
import argparse
import glob
import numpy as np
import shutil
import subprocess
import multiprocessing as mp
from functools import partial
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from skimage.io import imsave, imread
import pyvista as pv

from cc3d import CompuCellSetup
from cc3d.core.XMLUtils import ElementCC3D
from cc3d.core.PySteppables import SteppableBasePy


# ---------------------------
# Constants
# ---------------------------
FORCE_AMPLITUDE = 6.5  # Fixed amplitude for all experiments
NUM_RUNS = 100  # Number of runs per condition for statistics
NUM_WORKERS = 10 # min(mp.cpu_count(), NUM_RUNS)  # Number of parallel workers
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
    - Records trajectory for MSD calculation
    """

    def __init__(
        self,
        frequency=1,
        force_amplitude=FORCE_AMPLITUDE,
        direction_change_interval=DIRECTION_CHANGE_INTERVAL,
        outputdir="./output/",
        record_frequency=10,
        snapshot_frequency=SNAPSHOT_FREQUENCY
    ):
        super().__init__(frequency)

        self.force_amplitude = force_amplitude
        self.direction_change_interval = direction_change_interval
        self.outputdir = outputdir
        self.record_frequency = record_frequency
        self.snapshot_frequency = snapshot_frequency

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

    def save_instance_segmentation(self, mcs):
        """Save the current cell field state as a TIFF image."""
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

        # Save snapshot at regular intervals
        if mcs % self.snapshot_frequency == 0:
            self.save_instance_segmentation(mcs)

    def finish(self):
        df = pd.DataFrame(self.trajectories)
        trajectory_file = os.path.join(self.outputdir, "trajectories.csv")
        df.to_csv(trajectory_file, index=False)
        print(f"Trajectories saved to {trajectory_file}")


# ---------------------------
# Analysis Functions
# ---------------------------

def calculate_msd(trajectory_file):
    """
    Calculate Mean Square Displacement from trajectory data.
    MSD(t) = <|r(t) - r(0)|^2>

    Returns:
        tuple: (times, msd_values)
    """
    df = pd.read_csv(trajectory_file)

    if len(df) < 2:
        return np.array([0]), np.array([0])

    # Get initial position
    x0 = df.iloc[0]['x']
    y0 = df.iloc[0]['y']
    z0 = df.iloc[0]['z']

    # Calculate squared displacement for each time point
    df['sq_disp'] = (df['x'] - x0)**2 + (df['y'] - y0)**2 + (df['z'] - z0)**2

    # Time (MCS) and MSD
    times = df['mcs'].values
    msd = df['sq_disp'].values

    return times, msd


def power_law(t, D, alpha):
    """Power law function: MSD = D * t^alpha"""
    return D * np.power(t, alpha)


def fit_power_law(times, msd):
    """
    Fit MSD data to power law: MSD = D * t^alpha
    
    For normal diffusion: alpha = 1
    For subdiffusion (confined): alpha < 1
    For superdiffusion (directed): alpha > 1
    
    Returns:
        tuple: (D, alpha, D_err, alpha_err, r_squared)
    """
    if len(times) < 3:
        return 0.0, 1.0, 0.0, 0.0, 0.0
    
    # Remove t=0 point to avoid log(0) issues
    mask = times > 0
    t_fit = times[mask]
    msd_fit = msd[mask]
    
    if len(t_fit) < 2:
        return 0.0, 1.0, 0.0, 0.0, 0.0
    
    try:
        # Initial guess: D=1, alpha=1 (normal diffusion)
        popt, pcov = curve_fit(power_law, t_fit, msd_fit, p0=[1.0, 1.0], 
                               bounds=([0, 0], [np.inf, 3]), maxfev=5000)
        D, alpha = popt
        perr = np.sqrt(np.diag(pcov))
        D_err, alpha_err = perr
        
        # Calculate R²
        msd_pred = power_law(t_fit, D, alpha)
        ss_res = np.sum((msd_fit - msd_pred)**2)
        ss_tot = np.sum((msd_fit - np.mean(msd_fit))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return D, alpha, D_err, alpha_err, r_squared
    except Exception as e:
        print(f"Power law fit failed: {e}")
        return 0.0, 1.0, 0.0, 0.0, 0.0


def load_trajectory(trajectory_file):
    """Load trajectory data from CSV file."""
    df = pd.read_csv(trajectory_file)
    return df


def generate_3d_movie(run_dir, output_filename="cell_movie.gif", fps=10, lattice_size=LATTICE_SIZE):
    """
    Generate a 3D movie from TIFF snapshots using PyVista.
    Shows the cell moving within the full lattice space.
    
    Args:
        run_dir: Directory containing snapshot_*.tiff files
        output_filename: Name of the output movie file (GIF format for compatibility)
        fps: Frames per second for the movie
        lattice_size: Size of the Potts lattice (for fixed camera view)
    """
    # Find all snapshot files
    snapshot_pattern = os.path.join(run_dir, "snapshot_*.tiff")
    snapshot_files = sorted(glob.glob(snapshot_pattern))
    
    if not snapshot_files:
        print(f"No snapshots found in {run_dir}")
        return None
    
    print(f"Generating 3D movie from {len(snapshot_files)} snapshots...")
    
    # Configure PyVista for off-screen rendering
    pv.OFF_SCREEN = True
    
    output_path = os.path.join(run_dir, output_filename)
    frames = []
    
    # Fixed camera position for consistent view of entire lattice
    # Camera at a distance to see the full lattice
    center = lattice_size / 2
    distance = lattice_size * 2.5
    camera_position = [
        (center + distance, center + distance, center + distance),  # Camera position
        (center, center, center),  # Focal point (center of lattice)
        (0, 0, 1)  # Up vector
    ]
    
    try:
        for idx, snapshot_file in enumerate(snapshot_files):
            try:
                # Read the TIFF stack (Z, Y, X)
                img = imread(snapshot_file)
                
                # Create binary mask (cell vs background)
                binary = (img > 0).astype(np.float32)
                
                if binary.sum() == 0:
                    print(f"  Skipping empty frame {idx}")
                    continue
                
                # Create plotter for this frame
                plotter = pv.Plotter(off_screen=True, window_size=(800, 800))
                
                # Add wireframe box showing lattice boundaries
                box = pv.Box(bounds=(0, lattice_size, 0, lattice_size, 0, lattice_size))
                plotter.add_mesh(box, style='wireframe', color='gray', 
                                line_width=1, opacity=0.3)
                
                # Add axes at origin for reference
                plotter.add_axes(line_width=2, labels_off=False)
                
                # Transpose to (X, Y, Z) for PyVista and create uniform grid
                binary_xyz = np.transpose(binary, (2, 1, 0))  # (Z,Y,X) -> (X,Y,Z)
                
                # Create ImageData with point data
                grid = pv.ImageData(dimensions=np.array(binary_xyz.shape))
                grid.point_data["values"] = binary_xyz.flatten(order="F")
                
                # Extract the cell surface using marching cubes via contour
                contour = grid.contour([0.5], scalars="values")
                
                if contour.n_points > 0:
                    plotter.add_mesh(contour, color='dodgerblue', 
                                    smooth_shading=True, opacity=0.9)
                else:
                    # If contour fails, just show the occupied voxels as points
                    points = np.argwhere(binary_xyz > 0)
                    if len(points) > 0:
                        cloud = pv.PolyData(points)
                        plotter.add_mesh(cloud, color='dodgerblue', 
                                        point_size=5, render_points_as_spheres=True)
                
                # Set fixed camera and background
                plotter.set_background('white')
                plotter.camera_position = camera_position
                
                # Add frame number text
                mcs = int(os.path.basename(snapshot_file).replace('snapshot_', '').replace('.tiff', ''))
                plotter.add_text(f"MCS: {mcs}", position='upper_left', font_size=14)
                
                # Capture frame as image
                frame = plotter.screenshot(return_img=True)
                frames.append(frame)
                plotter.close()
                
                if (idx + 1) % 10 == 0:
                    print(f"  Processed {idx + 1}/{len(snapshot_files)} frames")
                    
            except Exception as frame_error:
                print(f"  Error processing frame {idx}: {frame_error}")
                continue
        
        # Save as GIF using imageio
        if frames:
            import imageio
            imageio.mimsave(output_path, frames, fps=fps, loop=0)
            print(f"3D movie saved to {output_path}")
            return output_path
        else:
            print("No frames generated - check if TIFF files contain cell data")
            return None
        
    except Exception as e:
        print(f"Error generating 3D movie: {e}")
        import traceback
        traceback.print_exc()
        return None


# ---------------------------
# Multi-run plotting functions
# ---------------------------

def plot_combined_trajectories_3d(trajectory_files, output_dir, title_suffix="", lattice_size=LATTICE_SIZE):
    """Create 3D plot with all trajectories from multiple runs."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(trajectory_files)))
    
    for idx, traj_file in enumerate(trajectory_files):
        df = pd.read_csv(traj_file)
        ax.plot(df['x'], df['y'], df['z'], '-', linewidth=1, alpha=0.7,
                color=colors[idx], label=f'Run {idx}')
        ax.scatter(df.iloc[0]['x'], df.iloc[0]['y'], df.iloc[0]['z'],
                   c=[colors[idx]], s=50, marker='o')
        ax.scatter(df.iloc[-1]['x'], df.iloc[-1]['y'], df.iloc[-1]['z'],
                   c=[colors[idx]], s=50, marker='s')

    # Set axis limits to match Potts lattice dimensions
    ax.set_xlim(0, lattice_size)
    ax.set_ylim(0, lattice_size)
    ax.set_zlim(0, lattice_size)
    
    ax.set_xlabel('X (pixels)')
    ax.set_ylabel('Y (pixels)')
    ax.set_zlabel('Z (pixels)')
    ax.set_title(f'3D Cell Trajectories - {len(trajectory_files)} Runs{title_suffix}')
    ax.legend(loc='upper left', fontsize=8)

    plot_file = os.path.join(output_dir, "trajectories_3d_combined.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Combined 3D trajectory plot saved to {plot_file}")


def plot_combined_msd(trajectory_files, output_dir, title_suffix=""):
    """
    Create MSD plot with all runs, mean, and std deviation.
    Fits to power law MSD = D * t^alpha to characterize diffusion type.
    """
    all_times = []
    all_msd = []
    
    # Collect MSD data from all runs
    for traj_file in trajectory_files:
        times, msd = calculate_msd(traj_file)
        all_times.append(times)
        all_msd.append(msd)
    
    # Find common time points (assuming all runs have same time points)
    times = all_times[0]
    msd_matrix = np.array(all_msd)
    
    # Calculate mean and std
    mean_msd = np.mean(msd_matrix, axis=0)
    std_msd = np.std(msd_matrix, axis=0)
    
    # Fit power law to mean MSD: MSD = D * t^alpha
    D_eff, alpha_eff, D_err, alpha_err, r_squared = fit_power_law(times, mean_msd)
    
    # Calculate D and alpha for each run
    D_values = []
    alpha_values = []
    for msd in all_msd:
        D, alpha, _, _, _ = fit_power_law(times, msd)
        D_values.append(D)
        alpha_values.append(alpha)
    D_std = np.std(D_values)
    D_mean = np.mean(D_values)
    alpha_std = np.std(alpha_values)
    alpha_mean = np.mean(alpha_values)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(trajectory_files)))
    
    # Plot individual runs
    for idx, msd in enumerate(all_msd):
        ax.plot(times, msd, '-', linewidth=1, alpha=0.4, color=colors[idx], label=f'Run {idx}')
    
    # Plot mean with error band
    ax.plot(times, mean_msd, 'k-', linewidth=2.5, label='Mean')
    ax.fill_between(times, mean_msd - std_msd, mean_msd + std_msd, 
                    alpha=0.3, color='gray', label='±1 Std Dev')
    
    # Add power law fit line
    t_fit = times[times > 0]
    fit_line = power_law(t_fit, D_eff, alpha_eff)
    ax.plot(t_fit, fit_line, 'r--', linewidth=2, alpha=0.7, 
            label=f'Power law fit: D·t^α\nD={D_mean:.4f}±{D_std:.4f}\nα={alpha_mean:.3f}±{alpha_std:.3f}')
    
    ax.set_xlabel('Time (MCS)', fontsize=14)
    ax.set_ylabel('MSD (pixels²)', fontsize=14)
    ax.set_title(f'Mean Square Displacement - {len(trajectory_files)} Runs{title_suffix}\n'
                 f'D = {D_mean:.4f} ± {D_std:.4f}, α = {alpha_mean:.3f} ± {alpha_std:.3f} (R² = {r_squared:.4f})', 
                 fontsize=14)
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = os.path.join(output_dir, "msd_combined.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Combined MSD plot saved to {plot_file}")
    
    return D_mean, D_std, alpha_mean, alpha_std, r_squared


# ---------------------------
# Simulation Runners
# ---------------------------

def run_single_cc3d_simulation(
    elastic_modulus=0.0,
    relaxation_time=0.0,
    output_dir="./output",
    run_index=0
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
            record_frequency=10
        )
    )

    # Run simulation
    CompuCellSetup.run()

    return run_dir


def _run_single_worker(args):
    """
    Worker function for parallel execution of single runs.
    Args is a tuple: (run_idx, elastic_modulus, relaxation_time, output_dir, script_path, python_exe)
    """
    run_idx, elastic_modulus, relaxation_time, output_dir, script_path, python_exe = args
    
    print(f"  Starting run {run_idx} (elastic={elastic_modulus}, tau={relaxation_time})")
    
    # Launch simulation in subprocess
    cmd = [
        python_exe, script_path, '--run_single',
        '--elastic', str(elastic_modulus),
        '--tau', str(relaxation_time),
        '--output', output_dir,
        '--run_index', str(run_idx)
    ]
    
    result = subprocess.run(cmd, capture_output=True)
    
    if result.returncode != 0:
        print(f"  WARNING: Run {run_idx} failed!")
        return None
    
    run_dir = os.path.join(output_dir, f"run_{run_idx}")
    traj_file = os.path.join(run_dir, "trajectories.csv")
    
    if os.path.exists(traj_file):
        print(f"  Completed run {run_idx}")
        # Generate 3D movie for this run
        # generate_3d_movie(run_dir)
        return traj_file
    
    return None


def run_multi_simulation(
    elastic_modulus=0.0,
    relaxation_time=0.0,
    output_base_dir="./single_test",
    num_runs=NUM_RUNS,
    num_workers=NUM_WORKERS
):
    """
    Run multiple simulations for statistical analysis.
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
        (run_idx, elastic_modulus, relaxation_time, output_dir, script_path, python_exe)
        for run_idx in range(num_runs)
    ]

    # Run in parallel using multiprocessing Pool
    with mp.Pool(processes=num_workers) as pool:
        results = pool.map(_run_single_worker, worker_args)

    # Collect successful trajectory files
    trajectory_files = [f for f in results if f is not None]
    print(f"\nCompleted {len(trajectory_files)}/{num_runs} runs successfully")

    # Generate combined plots
    if trajectory_files:
        title_suffix = f"\n(elastic={elastic_modulus}, tau={relaxation_time})"
        plot_combined_trajectories_3d(trajectory_files, output_dir, title_suffix)
        D_mean, D_std, alpha_mean, alpha_std, r_squared = plot_combined_msd(trajectory_files, output_dir, title_suffix)
        
        # Save summary
        summary = {
            'elastic_modulus': elastic_modulus,
            'relaxation_time': relaxation_time,
            'num_runs': len(trajectory_files),
            'D_eff_mean': D_mean,
            'D_eff_std': D_std,
            'alpha_mean': alpha_mean,
            'alpha_std': alpha_std,
            'r_squared': r_squared
        }
        summary_df = pd.DataFrame([summary])
        summary_df.to_csv(os.path.join(output_dir, "diffusion_summary.csv"), index=False)
        
        print(f"\n{'='*60}")
        print(f"RESULTS: D = {D_mean:.6f} ± {D_std:.6f} pixels²/MCS^α")
        print(f"         α = {alpha_mean:.4f} ± {alpha_std:.4f}")
        print(f"         R² = {r_squared:.4f}")
        print(f"{'='*60}\n")
        
        return output_dir, D_mean, D_std, alpha_mean, alpha_std
    
    return output_dir, 0.0, 0.0, 1.0, 0.0


def run_parameter_sweep():
    """
    Run a parameter sweep over elastic_modulus and relaxation_time.
    Each condition runs NUM_RUNS simulations for statistics.
    """
    output_base_dir = "./sweep_results"
    if os.path.exists(output_base_dir):
        shutil.rmtree(output_base_dir)
    os.makedirs(output_base_dir, exist_ok=True)

    # Parameters to sweep
    elastic_moduli = [0, 0.33, 0.66, 1.0, 1.33]
    relaxation_times = [0, 1000, 5000, 10000]

    results = []

    total_conditions = len(elastic_moduli) * len(relaxation_times)
    condition_count = 0

    for tau in relaxation_times:
        for elastic in elastic_moduli:
            condition_count += 1
            print(f"\n{'#'*60}")
            print(f"Condition {condition_count}/{total_conditions}")
            print(f"tau={tau}, elastic={elastic}")
            print(f"{'#'*60}\n")

            output_dir, D_mean, D_std, alpha_mean, alpha_std = run_multi_simulation(
                elastic_modulus=elastic,
                relaxation_time=tau,
                output_base_dir=output_base_dir,
                num_runs=NUM_RUNS
            )

            results.append({
                'tau': tau,
                'elastic_modulus': elastic,
                'D_eff_mean': D_mean,
                'D_eff_std': D_std,
                'alpha_mean': alpha_mean,
                'alpha_std': alpha_std
            })

    # Save results
    results_df = pd.DataFrame(results)
    results_file = os.path.join(output_base_dir, "diffusion_vs_parameters.csv")
    results_df.to_csv(results_file, index=False)
    print(f"\nResults saved to {results_file}")

    # Create summary plots
    create_summary_plots(results_df, output_base_dir)

    return results_df


def create_summary_plots(results_df, output_dir):
    """Create summary plots of diffusion coefficient and alpha vs parameters."""
    tau_values = sorted(results_df['tau'].unique())
    elastic_values = sorted(results_df['elastic_modulus'].unique())
    
    # Plot D vs elastic_modulus for different tau values
    colors = plt.cm.viridis(np.linspace(0, 1, len(tau_values)))

    plt.figure(figsize=(12, 7))
    for idx, tau in enumerate(tau_values):
        subset = results_df[results_df['tau'] == tau].sort_values('elastic_modulus')
        plt.errorbar(subset['elastic_modulus'], subset['D_eff_mean'],
                     yerr=subset['D_eff_std'], fmt='o-', linewidth=2, markersize=8,
                     color=colors[idx], markeredgecolor='black', capsize=5,
                     label=f'τ = {tau}')

    plt.xlabel('Elastic Modulus (λ_elastic)', fontsize=14)
    plt.ylabel('Diffusion Coefficient D (pixels²/MCS^α)', fontsize=14)
    plt.title('Diffusion Coefficient vs Elastic Modulus\n(MaxwellMedium)', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_file = os.path.join(output_dir, "D_vs_elastic_summary.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Summary plot saved to {plot_file}")

    # Plot alpha vs elastic_modulus for different tau values
    plt.figure(figsize=(12, 7))
    for idx, tau in enumerate(tau_values):
        subset = results_df[results_df['tau'] == tau].sort_values('elastic_modulus')
        plt.errorbar(subset['elastic_modulus'], subset['alpha_mean'],
                     yerr=subset['alpha_std'], fmt='o-', linewidth=2, markersize=8,
                     color=colors[idx], markeredgecolor='black', capsize=5,
                     label=f'τ = {tau}')

    plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, label='Normal diffusion (α=1)')
    plt.xlabel('Elastic Modulus (λ_elastic)', fontsize=14)
    plt.ylabel('Anomalous Exponent α', fontsize=14)
    plt.title('Anomalous Exponent vs Elastic Modulus\n(MaxwellMedium: α<1 subdiffusion, α>1 superdiffusion)', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_file = os.path.join(output_dir, "alpha_vs_elastic_summary.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Summary plot saved to {plot_file}")

    # Plot D vs tau for different elastic values
    colors = plt.cm.plasma(np.linspace(0, 1, len(elastic_values)))

    plt.figure(figsize=(12, 7))
    for idx, elastic in enumerate(elastic_values):
        subset = results_df[results_df['elastic_modulus'] == elastic].sort_values('tau')
        plt.errorbar(subset['tau'], subset['D_eff_mean'],
                     yerr=subset['D_eff_std'], fmt='o-', linewidth=2, markersize=8,
                     color=colors[idx], markeredgecolor='black', capsize=5,
                     label=f'λ = {elastic}')

    plt.xlabel('Relaxation Time (τ)', fontsize=14)
    plt.ylabel('Diffusion Coefficient D (pixels²/MCS^α)', fontsize=14)
    plt.title('Diffusion Coefficient vs Relaxation Time\n(MaxwellMedium)', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_file = os.path.join(output_dir, "D_vs_tau_summary.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Summary plot saved to {plot_file}")

    # Plot alpha vs tau for different elastic values
    plt.figure(figsize=(12, 7))
    for idx, elastic in enumerate(elastic_values):
        subset = results_df[results_df['elastic_modulus'] == elastic].sort_values('tau')
        plt.errorbar(subset['tau'], subset['alpha_mean'],
                     yerr=subset['alpha_std'], fmt='o-', linewidth=2, markersize=8,
                     color=colors[idx], markeredgecolor='black', capsize=5,
                     label=f'λ = {elastic}')

    plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, label='Normal diffusion (α=1)')
    plt.xlabel('Relaxation Time (τ)', fontsize=14)
    plt.ylabel('Anomalous Exponent α', fontsize=14)
    plt.title('Anomalous Exponent vs Relaxation Time\n(MaxwellMedium: α<1 subdiffusion, α>1 superdiffusion)', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_file = os.path.join(output_dir, "alpha_vs_tau_summary.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Summary plot saved to {plot_file}")


# ---------------------------
# Main
# ---------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Test MaxwellMedium plugin with Brownian directed motion and diffusion analysis"
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
        default='./single_test',
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

    args = parser.parse_args()

    if args.sweep:
        print("Running parameter sweep...")
        print(f"Each condition will run {NUM_RUNS} simulations for statistics")
        results_df = run_parameter_sweep()
        print("\nParameter sweep completed!")
        print(results_df)
    elif args.run_single:
        # Internal: run single CC3D simulation (called by subprocess)
        run_single_cc3d_simulation(
            elastic_modulus=args.elastic,
            relaxation_time=args.tau,
            output_dir=args.output,
            run_index=args.run_index
        )
    elif args.single:
        print(f"Running multi-simulation ({args.num_runs} runs):")
        print(f"  elastic_modulus = {args.elastic}")
        print(f"  tau = {args.tau}")
        print(f"  force_amplitude = {FORCE_AMPLITUDE} (fixed)")
        
        output_dir, D_mean, D_std, alpha_mean, alpha_std = run_multi_simulation(
            elastic_modulus=args.elastic,
            relaxation_time=args.tau,
            output_base_dir=args.output,
            num_runs=args.num_runs
        )
        
        print(f"\nSimulation completed! Output in {output_dir}")
        print(f"D = {D_mean:.6f} ± {D_std:.6f}, α = {alpha_mean:.4f} ± {alpha_std:.4f}")
    else:
        # Default: run multi-simulation test
        print(f"Running multi-simulation test ({NUM_RUNS} runs):")
        print(f"  elastic_modulus = {args.elastic}")
        print(f"  tau = {args.tau}")
        print(f"  force_amplitude = {FORCE_AMPLITUDE} (fixed)")
        
        output_dir, D_mean, D_std, alpha_mean, alpha_std = run_multi_simulation(
            elastic_modulus=args.elastic,
            relaxation_time=args.tau,
            output_base_dir=args.output,
            num_runs=NUM_RUNS
        )
        
        print(f"\nTest completed!")
        print(f"D = {D_mean:.6f} ± {D_std:.6f} pixels²/MCS^α")
        print(f"α = {alpha_mean:.4f} ± {alpha_std:.4f}")
