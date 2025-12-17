# 🐳 Container Build

This directory contains the files needed to build the `gdavid57/cc3d-compile` Docker image from source.

## Pre-built Image

The easiest way to get started is to pull the pre-built image from Docker Hub:

```bash
docker pull gdavid57/cc3d-compile
```

## Building from Source

If you need to modify the image (e.g., add custom plugins), you can build it locally.

### Prerequisites

- Docker 20.10+

### Build Instructions

From the **repository root**:

```bash
docker build \
    -t cc3d-compile \
    -f container_build/Dockerfile \
    --build-context DeveloperZone=cc3d_plugins \
    .
```

> **Note:** The `--build-context` flag maps the `cc3d_plugins/` directory to the `DeveloperZone` context expected by the Dockerfile.

### What the Image Contains

- **Ubuntu 22.04** base
- **Miniconda** with Python 3.10
- **CompuCell3D 4.6.0** compiled from source
- **Custom plugins**: BacteriaShapeBis, MaxwellMedium
- All scientific Python dependencies (numpy, scipy, pandas, scikit-image, etc.)

## Running the Container

### Interactive Mode

```bash
docker run -it --rm gdavid57/cc3d-compile
```

### With Volume Mounts (for running benchmarks)

```bash
docker run -it --rm \
    -v /path/to/benchmarks:/script_cc3d \
    -v /path/to/whole_cell_potts_model:/whole_cell_potts_model \
    gdavid57/cc3d-compile
```

### With GUI Support (Linux/X11)

To run graphical applications (Player5, Twedit5) from the container:

```bash
# Allow Docker to access X11
xhost +local:docker

# Run with display forwarding
docker run -it --rm \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    gdavid57/cc3d-compile
```

## Environment Details

| Variable | Value |
|----------|-------|
| Conda environment | `cc3d_compile` |
| Python version | 3.10 |
| Working directory | `/script_cc3d` |
| User | `cc3duser` |

