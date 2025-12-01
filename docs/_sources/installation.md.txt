# Installation

## Install via PyPI

To install **PYATB** from PyPI, ensure that `mpi4py` is available in your Python environment. We recommend installing it via conda:

```shell
conda install -c conda-forge mpi4py
pip install pyatb
```

## Install from Source

The latest development version of PYATB can be obtained from the official GitHub repository:

```shell
git clone https://github.com/pyatb/pyatb.git
```

### Prerequisites

Currently, PYATB is only supported on Linux systems and has not been tested on Windows or Mac systems. To use PYATB, the following prerequisites must be installed:

- Python 3.8 or newer
- C++ compiler
- [pybind11](https://pybind11.readthedocs.io/en/stable/index.html)
- [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [mpi4py](https://mpi4py.readthedocs.io/en/stable/mpi4py.html)
- [Matplotlib](https://matplotlib.org/)
- [ASE](https://ase-lib.org/)


### Building Wheel Files (.whl)

PYATB provides a convenient Docker-based workflow for building **manylinux-compliant wheel packages**.

The `build_wheels.sh` script will automatically create a Docker environment and compile the wheel files:

```shell
./build_wheels.sh
```

If you need to build wheels for a specific Python version, modify the `CIBW_BUILD` variable inside the script. For example:

```shell
CIBW_BUILD="cp313-manylinux_x86_64"
```

will generate wheel files for **Python 3.13**.
