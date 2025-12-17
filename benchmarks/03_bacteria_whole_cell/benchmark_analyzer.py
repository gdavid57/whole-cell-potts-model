import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path


def analyze_benchmarks(csv_path: Path):
    """
    Analyzes benchmark results from a CSV file to generate performance graphs.

    This function performs the following steps:
    1.  Loads the benchmark data from the specified CSV file.
    2.  Converts the 'execution_time' column from a time string to total seconds.
    3.  Identifies the default configuration by finding the mode (most frequent value) for each parameter.
    4.  Determines which parameters were benchmarked by checking for columns with multiple unique values.
    5.  For each benchmarked parameter, it generates and saves a plot showing its impact
        on execution time relative to the default configuration.

    Args:
        csv_path (str): The path to the input CSV file.
    """
    # --- 1. Load Data and Pre-process ---
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: The file '{csv_path}' was not found.")
        return

    # Convert execution time to a numerical format (total seconds)
    # The 'errors='coerce'' will turn any parsing errors into NaT (Not a Time)
    df["execution_time_seconds"] = pd.to_timedelta(
        df["execution_time"], errors="coerce"
    ).dt.total_seconds()

    # Drop rows where time conversion failed, if any
    df.dropna(subset=["execution_time_seconds"], inplace=True)

    print("Data loaded and execution time converted to seconds.")
    print(f"Total number of runs found: {len(df)}")

    # --- 2. Find the Default Configuration ---
    # Identify parameter columns (exclude output/metadata columns)
    # We consider all columns except the time-related ones as potential parameters.
    # You can add other non-parameter columns to this list if needed.
    non_param_cols = [
        "execution_time",
        "execution_time_seconds",
        "time_reference",
        "outputdir",
        "trna_data_path",
        "mrna_data_path",
    ]
    param_cols = [col for col in df.columns if col not in non_param_cols]

    # The default configuration is the one with the most common value for each parameter.
    # We compute the mode for all parameter columns.
    default_config = df[param_cols].mode().iloc[0]

    print("\n--- Default Configuration ---")
    print(default_config)
    print("-----------------------------\n")

    # --- 3. Identify Benchmark Parameters ---
    # Find parameters that have more than one unique value
    varying_params = []
    for col in param_cols:
        if df[col].nunique() > 1:
            varying_params.append(col)

    if not varying_params:
        print(
            "No varying parameters found. All runs seem to have the same configuration."
        )
        return

    print(f"Found varying parameters to plot: {varying_params}\n")

    # Create a directory to save plots
    output_dir = f"benchmark_plots_{csv_path.stem}"
    os.makedirs(output_dir, exist_ok=True)
    print(f"Plots will be saved in the '{output_dir}' directory.")

    # --- 4. Generate a Plot for Each Varying Parameter ---
    sns.set_theme(style="whitegrid")

    for param in varying_params:
        # For a given parameter, we want to select all rows where the *other* parameters
        # match the default configuration. This isolates the effect of the current parameter.
        other_params = default_config.index.drop(param).tolist()

        # Create a boolean mask for rows that match the default config on all other parameters
        mask = (df[other_params] == default_config[other_params]).all(axis=1)

        # Filter the dataframe to get the data for our plot
        plot_data = df[mask].copy()

        # Check if we have enough data to plot
        if plot_data.empty or len(plot_data[param].unique()) < 2:
            print(
                f"--> Skipping plot for '{param}': Not enough data points where other params are default."
            )
            continue

        # Sort values for a clean plot
        plot_data.sort_values(by=param, inplace=True)

        plt.figure(figsize=(10, 6))

        # Plot the data points
        ax = sns.lineplot(
            data=plot_data,
            x=param,
            y="execution_time_seconds",
            marker="o",
            errorbar=None,
        )

        # Add labels for each data point
        for index, row in plot_data.iterrows():
            ax.text(
                row[param],
                row["execution_time_seconds"],
                f" {row['execution_time_seconds']:.0f}s",
                ha="left",
                va="center",
                fontsize=9,
            )

        plt.title(f'Impact of "{param}" on Execution Time', fontsize=16)
        plt.xlabel(param, fontsize=12)
        plt.ylabel("Execution Time (seconds)", fontsize=12)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()  # Adjust layout to make room for labels

        # Save the figure
        plot_filename = os.path.join(output_dir, f"benchmark_{param}.png")
        plt.savefig(plot_filename)
        plt.close()

        print(f"--> Generated plot for '{param}' and saved as '{plot_filename}'")


if __name__ == "__main__":
    # Assuming CSV files are in the same directory as the script.
    for filename in [
        "results_multithreading_wc_script.csv",
        "results_multiprocessing_wc_script.csv",
        "results_automated_wc_cc3d_script.csv",
    ]:
        path = Path(filename)

        if path.exists():
            analyze_benchmarks(path)
        else:
            print(f"{path} doesn't exist")
