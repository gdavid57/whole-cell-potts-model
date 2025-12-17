import os
import sys
import argparse
import random
import numpy as np
import json
from skimage.io import imsave

from cc3d import CompuCellSetup
from cc3d.core.XMLUtils import ElementCC3D
from cc3d.core.PySteppables import SteppableBasePy, MitosisSteppableBase
from cc3d.cpp import CompuCellExtraModules

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
_WHOLE_CELL_MODEL_DIR = os.path.join(_REPO_ROOT, "whole_cell_model")
if _WHOLE_CELL_MODEL_DIR not in sys.path:
    sys.path.insert(0, _WHOLE_CELL_MODEL_DIR)

import py_cell_state


# ---------------------------
# XML-like configuration
# ---------------------------
def configure_simulation(
    number_of_processors=1,
    debug_output_frequency=100000,
    potts_dimensions={"x": 128, "y": 128, "z": 128},
    potts_steps=100000,
    potts_temperature=37.0,
    neighbor_order=1,
    bacteria_rotation_criterion=1000,
    bacteria_random_amplitude=0.0,
    bacteria_lambda_shape=100.0,
):
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

    # Plugins
    cc3d.ElementCC3D("Plugin", {"Name": "MomentOfInertia"})
    cc3d.ElementCC3D("Plugin", {"Name": "Volume"})
    celltype = cc3d.ElementCC3D("Plugin", {"Name": "CellType"})
    celltype.ElementCC3D("CellType", {"TypeId": 0, "TypeName": "Medium"})
    celltype.ElementCC3D("CellType", {"TypeId": 1, "TypeName": "ECOLI"})
    cc3d.ElementCC3D("Plugin", {"Name": "CenterOfMass"})
    cc3d.ElementCC3D("Plugin", {"Name": "NeighborTracker"})

    # Plugin BacteriaShapeBis
    bact = cc3d.ElementCC3D("Plugin", {"Name": "BacteriaShapeBis"})
    bact.ElementCC3D("RotationCriterion", {}, bacteria_rotation_criterion)
    bact.ElementCC3D("RandomAmplitude", {}, bacteria_random_amplitude)
    bact.ElementCC3D("LambdaShape", {}, bacteria_lambda_shape)

    CompuCellSetup.setSimulationXMLDescription(cc3d)


# ---------------------------
# Steppable
# ---------------------------


class InitialLayoutSteppable(SteppableBasePy):
    def __init__(self, frequency=1, n_cells=1, radius=3, margin=15):
        super().__init__(frequency)

        self.n_cells = n_cells
        self.radius = radius
        self.margin = margin

    def start(self):
        for _ in range(self.n_cells):
            cx = random.randint(self.margin, self.dim.x - self.margin - 1)
            cy = random.randint(self.margin, self.dim.y - self.margin - 1)
            cz = random.randint(self.margin, self.dim.z - self.margin - 1)

            cell = self.new_cell(self.ECOLI)

            for dx in range(-self.radius, self.radius + 1):
                for dy in range(-self.radius, self.radius + 1):
                    for dz in range(-self.radius, self.radius + 1):
                        x = max(0, min(self.dim.x - 1, cx + dx))
                        y = max(0, min(self.dim.y - 1, cy + dy))
                        z = max(0, min(self.dim.z - 1, cz + dz))
                        self.cell_field[x, y, z] = cell


class GrowthSteppable(SteppableBasePy):
    def __init__(
        self,
        frequency=1,
        initial_major_length=15.0,
        minor_length=10.0,
        lambda_volume=100.0,
        cc3d_steps=10,
        outputdir="./snapshots/",
        snapshot_steps=None,
        trna_data_path="data/",
        mrna_data_path="data/LB/",
        init_cell_volume=3.5e-18,
        num_ribosomes=50000,
        time_reference=1200,
        whole_cell_steps=10000,
        step_to_stationary_stage=500000,
    ):
        super().__init__(frequency)

        # CC3D parameters
        self.initial_major_length = initial_major_length
        self.minor_length = minor_length
        self.lambda_volume = lambda_volume
        self.cc3d_steps = cc3d_steps
        self.outputdir = outputdir
        self.snapshot_steps = snapshot_steps
        os.makedirs(self.outputdir, exist_ok=True)

        # Whole Cell parameters
        self.trna_data_path = trna_data_path
        self.mrna_data_path = mrna_data_path
        self.init_cell_volume = init_cell_volume
        self.num_ribosomes = num_ribosomes
        self.time_reference = time_reference
        self.whole_cell_steps = whole_cell_steps
        self.step_to_stationary_stage = step_to_stationary_stage

        # Useful quantities
        self.bacteria_model_plugin = None
        self.growth_scale = self.initial_major_length / self.time_reference
        self.cc3d_initial_volume = self.spherocylinder_volume(
            self.initial_major_length, self.minor_length
        )
        self.volume_conversion_factor = self.init_cell_volume / self.cc3d_initial_volume
        self.initial_n_cells = 0

    def growth_model(self, lifetime):
        return self.growth_scale * lifetime

    def volume_conversion(self, cell_volume):
        return self.volume_conversion_factor * cell_volume

    def spherocylinder_volume(self, major_length, minor_length):
        return (
            0.25 * np.pi * minor_length**2 * (major_length - minor_length)
            + np.pi * minor_length**3 / 6.0
        )

    def save_instance_segmentation(self, mcs):
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

        imsave(self.outputdir + name, img.astype(np.uint8), check_contrast=False)

    def initialize_cell(self, cell):
        cell.lambdaVolume = self.lambda_volume
        cell.targetVolume = self.cc3d_initial_volume

        self.bacteria_model_plugin.setMajorAxisLength(cell, self.initial_major_length)
        self.bacteria_model_plugin.setMinorAxisLength(cell, self.minor_length)

        cell.dict["cell_state"] = py_cell_state.CellState(
            self.trna_data_path,
            self.mrna_data_path,
            self.init_cell_volume,
            self.num_ribosomes,
        )

        initial_state = cell.dict["cell_state"].run_step(
            self.initial_major_length, self.step_to_stationary_stage
        )

        cell.dict["producted_proteins"] = 0
        cell.dict["lifetime"] = 0
        cell.dict["lifetime_offset"] = initial_state.physical_time
        cell.dict["has_divided"] = False

    def start(self):
        self.bacteria_model_plugin = CompuCellExtraModules.getBacteriaShapeBisPlugin()

        for cell in self.cell_list:
            self.initialize_cell(cell)

        self.initial_n_cells = len(self.cell_list)

    def step(self, mcs: int):
        # FOR TESTING PURPOSE
        if len(self.cell_list) >= 4 * self.initial_n_cells:
            sys.exit(1)

        print(f"{len(self.cell_list)=}, {mcs=}, {self.cc3d_steps=}")

        for cell in self.cell_list:
            if cell.dict["has_divided"]:
                self.initialize_cell(cell)
                cell.dict["has_divided"] = False

            if mcs % self.cc3d_steps == 0:
                converted_volume = self.volume_conversion(cell.volume)

                result = cell.dict["cell_state"].run_step(
                    converted_volume, self.whole_cell_steps
                )
                cell.dict["producted_proteins"] += result.produced_proteins
                cell.dict["lifetime"] = result.physical_time

                new_major_length = self.initial_major_length + self.growth_model(
                    cell.dict["lifetime"] - cell.dict["lifetime_offset"]
                )

                self.bacteria_model_plugin.setMajorAxisLength(cell, new_major_length)
                cell.targetVolume = self.spherocylinder_volume(
                    new_major_length, self.minor_length
                )

                if new_major_length > 2 * self.initial_major_length:
                    self.shared_steppable_vars["cells_to_divide"].append(cell)

                print(
                    "MCS:",
                    mcs,
                    " Cell ID:",
                    cell.id,
                    " Volume:",
                    cell.volume,
                    " Target:",
                    cell.targetVolume,
                    " Major Length:",
                    new_major_length,
                    " COM:",
                    (cell.xCOM, cell.yCOM, cell.zCOM),
                    " Proteins:",
                    cell.dict["producted_proteins"],
                    " Lifetime:",
                    cell.dict["lifetime"],
                )

        if self.snapshot_steps and (mcs % self.snapshot_steps == 0):
            self.save_instance_segmentation(mcs)


class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        super().__init__(frequency)

    def start(self):
        self.shared_steppable_vars["cells_to_divide"] = []

    def step(self, mcs: int):
        for cell in self.shared_steppable_vars["cells_to_divide"]:
            self.divide_cell_along_minor_axis(cell)

        self.shared_steppable_vars["cells_to_divide"] = []

    def update_attributes(self):
        self.parent_cell.dict["has_divided"] = True
        self.child_cell.dict["has_divided"] = True


def run_simulation(config: dict):
    """Configuration & Run.

    Args:
        config (dict): configuration from JSON file.
    """
    # Potts & Bacteria configuration
    configure_simulation(
        number_of_processors=1, debug_output_frequency=10, **config["simulation"]
    )

    # Lattice Initialization Steppable
    CompuCellSetup.register_steppable(
        steppable=InitialLayoutSteppable(
            frequency=1, **config["lattice_initialization_steppable"]
        )
    )

    # Growth Steppable
    CompuCellSetup.register_steppable(
        steppable=GrowthSteppable(frequency=1, **config["growth_steppable"])
    )

    # Mitosis Steppable
    CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

    # Simulation Run
    CompuCellSetup.run()


def parse_json(input_string: str):
    """Parses JSON from a file or a string."""
    if os.path.exists(input_string):
        # It's a file path
        try:
            with open(input_string, "r") as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            print(
                f"Error: Could not decode JSON from file '{input_string}'. Details: {e}"
            )
            sys.exit(1)
    else:
        # It's a JSON string
        try:
            return json.loads(input_string)
        except json.JSONDecodeError as e:
            print(f"Error: Could not decode JSON from string. Details: {e}")
            sys.exit(1)


# ---------------------------
# Main: Configuration & Run
# ---------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse JSON from a file path or a JSON string."
    )
    parser.add_argument(
        "input",
        type=str,
        help="A file path to a JSON file or a JSON string.",
    )

    args = parser.parse_args()
    config = parse_json(args.input)

    run_simulation(config)
