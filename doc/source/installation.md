# Installation

## Download

The latest version of the PYATB package can be obtained from the GitHub repository:

```shell
git clone https://github.com/pyatb/pyatb.git
```

## Prerequisites

Currently, PYATB is only supported on Linux systems and has not been tested on Windows or Mac systems. To use PYATB, the following prerequisites must be installed:

- Python 3.7 or newer
- NumPy
- SciPy
- mpi4py
- Matplotlib
- C++ compiler
- Intel MKL 

## Install

You can install PYATB with the following command:

```shell
python setup.py install --record log
```

To customize the `setup.py` file, you must make changes to the **CXX** and **LAPACK_DIR** variables in line with your environment. **CXX** denotes the C++ compiler you intend to use, for instance, icpc (note that it should not be the mpi version). Furthermore, **LAPACK_DIR** is used to specify the Intel MKL path.

After completing the installation process, you can access the `pyatb` executable and corresponding module, which can be imported using the `import pyatb` command.

To uninstall the software, you simply delete the files generated during the Python installation process. A straightforward way to do this is to execute the following uninstallation command:

```shell
cat log | xargs rm -rf
```

**NOTE**: Be sure to delete all files related to pyatb. Otherwise it will cause unnecessary impact on the next installation.