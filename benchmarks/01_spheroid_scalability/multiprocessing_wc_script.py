from concurrent.futures import Future, ProcessPoolExecutor
import os
import sys
import argparse
import numpy as np
import json
import random
import time
from skimage.io import imsave

from cc3d import CompuCellSetup
from cc3d.core.XMLUtils import ElementCC3D
from cc3d.core.PySteppables import SteppableBasePy, MitosisSteppableBase

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
    celltype.ElementCC3D("CellType", {"TypeId": 1, "TypeName": "TUMOR"})
    cc3d.ElementCC3D("Plugin", {"Name": "CenterOfMass"})
    cc3d.ElementCC3D("Plugin", {"Name": "NeighborTracker"})

    contact = cc3d.ElementCC3D("Plugin", {"Name": "Contact"})
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Medium"}, 0.0)
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "TUMOR"}, 10.0)
    contact.ElementCC3D("Energy", {"Type1": "TUMOR", "Type2": "TUMOR"}, 5.0)
    contact.ElementCC3D("NeighborOrder", {}, neighbor_order)

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
        cx = self.dim.x // 2
        cy = self.dim.y // 2
        cz = self.dim.z // 2

        for _ in range(self.n_cells):
            cell = self.new_cell(self.TUMOR)

            for dx in range(-self.radius, self.radius + 1):
                for dy in range(-self.radius, self.radius + 1):
                    for dz in range(-self.radius, self.radius + 1):
                        x = max(0, min(self.dim.x - 1, cx + dx))
                        y = max(0, min(self.dim.y - 1, cy + dy))
                        z = max(0, min(self.dim.z - 1, cz + dz))
                        self.cell_field[x, y, z] = cell


def run_cell_state_step(cell_state, volume: float, steps: int):
    """
    Function wrapper to be used in ProcessPoolExecutor  to be able to serialise
    the new state of the CellState object.
    """
    cell_state_step_result = cell_state.run_step(volume, steps)

    return cell_state, cell_state_step_result


class GrowthSteppableThreaded(SteppableBasePy):
    def __init__(
        self,
        frequency=1,
        lambda_volume=100.0,
        target_volume_min=500.0,
        target_volume_max=1000.0,
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
        generation_csv_path=None,
    ):
        super().__init__(frequency)

        # CC3D parameters
        self.lambda_volume = lambda_volume
        self.target_volume_min = target_volume_min
        self.target_volume_max = target_volume_max
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
        self.volume_growth_rate = (
            self.target_volume_max - self.target_volume_min
        ) / self.time_reference
        self.cc3d_initial_volume = self.target_volume_min
        self.volume_conversion_factor = self.init_cell_volume / self.cc3d_initial_volume
        self.initial_n_cells = 0
        self.generation_csv_path = (
            generation_csv_path
            if generation_csv_path is not None
            else os.path.join(self.outputdir, "generation_times.csv")
        )
        self._start_time = None
        self._next_doubling_target = None

    def growth_model(self, lifetime):
        return min(
            self.target_volume_max,
            self.target_volume_min + self.volume_growth_rate * lifetime,
        )

    def volume_conversion(self, cell_volume):
        return self.volume_conversion_factor * cell_volume

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
        imsave(self.outputdir + name, img.astype(np.uint16), check_contrast=False)

    def initialize_divided_cells(self):
        divided_cells = [cell for cell in self.cell_list if cell.dict["has_divided"]]

        with ProcessPoolExecutor() as executor:
            futures: list[Future] = []

            for cell in divided_cells:
                cell.lambdaVolume = self.lambda_volume
                cell.targetVolume = self.target_volume_min

                cell_state = py_cell_state.CellState(
                    self.trna_data_path,
                    self.mrna_data_path,
                    self.init_cell_volume,
                    self.num_ribosomes,
                )

                cell.dict["cell_state"] = cell_state

                future = executor.submit(
                    run_cell_state_step,
                    cell_state,
                    self.volume_conversion(self.target_volume_min),
                    self.step_to_stationary_stage,
                )
                futures.append(future)

            for future, cell in zip(futures, divided_cells):
                new_cell_state, initial_state = future.result()
                cell.dict["cell_state"] = new_cell_state
                cell.dict["producted_proteins"] = 0
                cell.dict["lifetime"] = 0
                cell.dict["lifetime_offset"] = initial_state.physical_time
                cell.dict["has_divided"] = False
                cell.dict["target_volume"] = self.target_volume_min

    def start(self):
        for cell in self.cell_list:
            cell.dict["has_divided"] = True

        self.initial_n_cells = len(self.cell_list)
        self._start_time = time.time()
        self._next_doubling_target = max(2, self.initial_n_cells * 2)

        # reset CSV
        with open(self.generation_csv_path, "w", encoding="utf-8") as csv_file:
            csv_file.write("cell_count,elapsed_seconds,mcs\n")

    def process_cells(self):
        results = []

        with ProcessPoolExecutor() as executor:
            futures: list[Future] = []

            for cell in self.cell_list:
                converted_volume = self.volume_conversion(cell.volume)

                # Exécuter le modèle Whole Cell
                future = executor.submit(
                    run_cell_state_step,
                    cell.dict["cell_state"],
                    converted_volume,
                    self.whole_cell_steps,
                )
                futures.append(future)

            for future, cell in zip(futures, self.cell_list):
                new_cell_state, result = future.result()
                cell.dict["cell_state"] = new_cell_state
                cell.dict["producted_proteins"] += result.produced_proteins
                cell.dict["lifetime"] = result.physical_time

                lifetime_in_cycle = cell.dict["lifetime"] - cell.dict["lifetime_offset"]
                new_target_volume = self.growth_model(lifetime_in_cycle)
                should_divide = cell.volume >= self.target_volume_max

                results.append(
                    {
                        "cell": cell,
                        "new_target_volume": new_target_volume,
                        "should_divide": should_divide,
                        "success": True,
                    }
                )
                # except Exception as e:
                #     print(f"Error processing cell {cell.id}: {e}")
                #     results.append({"cell": cell, "success": False, "error": str(e)})

        return results

    def _record_generation_if_needed(self, mcs: int):
        if self._start_time is None or self._next_doubling_target is None:
            return

        while len(self.cell_list) >= self._next_doubling_target:
            elapsed = time.time() - self._start_time
            with open(self.generation_csv_path, "a", encoding="utf-8") as csv_file:
                csv_file.write(
                    f"{len(self.cell_list)},{elapsed:.6f},{mcs}\n"
                )

            self._next_doubling_target *= 2

    def step(self, mcs: int):

        self._record_generation_if_needed(mcs)

        # restore legacy stop condition used previously (testing safeguard)
        if len(self.cell_list) >= 1024:
            sys.exit(1)

        self.initialize_divided_cells()

        if mcs % self.cc3d_steps == 0:
            cells_to_divide = []

            results = self.process_cells()

            for result in results:
                if result["success"]:
                    cell = result["cell"]

                    cell.targetVolume = result["new_target_volume"]
                    cell.dict["target_volume"] = result["new_target_volume"]

                    if result["should_divide"]:
                        cells_to_divide.append(cell)

                else:
                    print(
                        f"Failed to process cell: {result.get('error', 'Unknown error')}"
                    )

            if cells_to_divide:
                self.shared_steppable_vars["cells_to_divide"] = cells_to_divide

        if self.snapshot_steps and (mcs % self.snapshot_steps == 0):
            self.save_instance_segmentation(mcs)

        for cell in self.cell_list:
            if cell.id % 50 == 1:
                print(
                    "MCS:",
                    mcs,
                    " Cell ID:",
                    cell.id,
                    " Volume:",
                    cell.volume,
                    " Target:",
                    cell.targetVolume,
                    " Target Dict:",
                    cell.dict.get("target_volume"),
                    " COM:",
                    (cell.xCOM, cell.yCOM, cell.zCOM),
                    " Proteins:",
                    cell.dict.get("producted_proteins"),
                    " Lifetime:",
                    cell.dict.get("lifetime"),
                )


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
        
        self.clone_parent_2_child()

        self.parent_cell.dict["has_divided"] = True
        self.child_cell.dict["has_divided"] = True


def run_simulation(config: dict):
    """Configuration & Run.

    Args:
        config (dict): configuration from JSON file.
    """
    # Potts configuration
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
        steppable=GrowthSteppableThreaded(frequency=1, **config["growth_steppable"])
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
