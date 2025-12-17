import csv
import json
import subprocess
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

PYTHON_SCRIPT_PATH_LIST = [
    Path("multiprocessing_wc_script.py"),
    Path("multithreading_wc_script.py"),
    Path("automated_wc_cc3d_script.py"),
]
DEFAULT_CONFIG_PATH = Path("default_config.json")
RANGES_TO_TEST_PATH = Path("ranges_to_test.json")

FIELD_NAMES = [
    "potts_dimension_x",
    "potts_dimension_y",
    "potts_dimension_z",
    "potts_steps",
    "potts_temperature",
    "neighbor_order",
    "bacteria_rotation_criterion",
    "bacteria_random_amplitude",
    "bacteria_lambda_shape",
    "n_cells",
    "radius",
    "margin",
    "initial_major_length",
    "minor_length",
    "lambda_volume",
    "cc3d_steps",
    "outputdir",
    "snapshot_steps",
    "trna_data_path",
    "mrna_data_path",
    "init_cell_volume",
    "num_ribosomes",
    "time_reference",
    "whole_cell_steps",
    "step_to_stationary_stage",
    "execution_time",
]


def time_simulation(python_script_path: Path, str_config: str) -> timedelta:
    """
    Run the simulation with subprocess to avoid the crash of this script
    when segfault error occurs in the executed script.
    This error happens at the very end of each run but does not prevent
    the data/output files generation. This error can be ignored for now.
    """
    start = datetime.now()

    try:
        subprocess.run(
            ["python", str(python_script_path), str_config],
            capture_output=True,
            text=True,
            check=True,
        )

    except Exception as exc:
        print(exc)

    finally:
        return datetime.now() - start


def flatten_config(config: dict) -> dict:
    potts_dimensions = {
        "potts_dimension_x": config["simulation"]["potts_dimensions"]["x"],
        "potts_dimension_y": config["simulation"]["potts_dimensions"]["y"],
        "potts_dimension_z": config["simulation"]["potts_dimensions"]["z"],
    }
    flatten_data = (
        potts_dimensions
        | config["simulation"]
        | config["lattice_initialization_steppable"]
        | config["growth_steppable"]
    )
    del flatten_data["potts_dimensions"]

    return flatten_data


def stringify(config: dict) -> str:
    return json.dumps(config, sort_keys=True)


def generate_str_config_set(
    default_config: dict[str, dict[str, Any]], ranges_to_test: dict[str, dict[str, Any]]
) -> set[str]:
    """
    Returns a set of string representations of the configuration dictionaries.
    They are derivated from the default configuration.
    """
    # Set to ensure unicity of each config
    str_config_set = set()

    # Add the default config to the set
    str_default_config = stringify(default_config)
    str_config_set.add(str_default_config)

    for key, subconfig in default_config.items():
        for parameter_name in subconfig:
            # Get the range to test for given parameter name
            parameter_to_test = ranges_to_test[key].get(parameter_name)

            if isinstance(parameter_to_test, list):
                for value in parameter_to_test:
                    # Create a deep copy from the default config
                    new_config = json.loads(str_default_config)
                    new_config[key][parameter_name] = value
                    str_config_set.add(stringify(new_config))

    return str_config_set


if __name__ == "__main__":
    with DEFAULT_CONFIG_PATH.open() as json_file:
        default_config: dict[str, dict] = json.load(json_file)

    with RANGES_TO_TEST_PATH.open() as json_file:
        ranges_to_test: dict[str, dict] = json.load(json_file)

    # Generate unique configurations
    str_config_set = generate_str_config_set(default_config, ranges_to_test)
    print(f"NUMBER OF GENERATED CONFIGURATIONS: {len(str_config_set)}")

    for script_path in PYTHON_SCRIPT_PATH_LIST:
        if script_path.exists():
            print(f"Run {script_path} script")
            output_filepath = Path(f"results_{script_path.stem}.csv")

            with output_filepath.open("w") as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=FIELD_NAMES)
                writer.writeheader()

            for str_config in str_config_set:
                time_delta = time_simulation(script_path, str_config)
                # Flatten the deserialized config to be writeable into the CSV file
                flattened_data = flatten_config(json.loads(str_config))
                # Add execution time
                flattened_data["execution_time"] = time_delta

                with output_filepath.open("a") as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=FIELD_NAMES)
                    writer.writerow(flattened_data)
