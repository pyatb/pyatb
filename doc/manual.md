
- [Capabilities](#capabilities)
- [Download and install](#download-and-install)
  - [Installation](#installation)
    - [Requirements](#requirements)
    - [Installation from source](#installation-from-source)
- [Quickstart guide](#quickstart-guide)
  - [Input files](#input-files)
- [BUG](#bug)


# Capabilities

PYATB provides the following functionalities:

1. band structue
2. band unfolding
3. fermi energy
4. fermi surface
5. find nodes
6. JDOS
7. PDOS
8. spin texture
9. AHC
10. berry curvature
11. berry curvature dipole
12. chern number
13. chirality
14. CPGE
15. optical conductivity
16. polarization
17. shift current
18. wilson loop

# Download and install

## Installation

### Requirements

- Python 3.7 or newer
- NumPy
- SciPy
- mpi4py
- Matplotlib
- C++ compiler
- BLAS, LAPACK

### Installation from source

Install the PYATB program into the python environment:

```shell {.line-numbers}
python setup.py install --record log
```

After the installation is successful, the executable program `pyatb` exists in the python environment.




# Quickstart guide

The PYATB program can be used as a general program, that is, read the specified `Input` file, and perform the calculation of each function described in it. After preparing the Input file and the corresponding HR, SR, and rR files, directly execute the `mpirun -np process_num pyatb` command to calculate.

## Input files

- The `Input` file

    This file describes the basic information of the system and the parameters required to calculate the function. For a complete list of the input parameters, please consult this [input list](input.md).

- The HR file

    This file describes the Hamiltonian of the system in real space R.

- The SR file

    This file describes the overlap matrix S in real space R.

- The rR file
  
    This file describes the dipole matrix r in real space R.

- The structure file

    This file describes the structural information about the system, e.g., lattice constant, lattice vectors, and positions of the atoms within a unit cell. The structure file format of ABACUS is currently used.

- The NAOs files

    This file describes the basis sets used to expand the Hamiltonian, the numerical atomic orbitals.



# BUG
1. "import matplotlib"在"from mpi4py import MPI"会报错，显示"/lib64/libz.so.1: version `ZLIB_1.2.9\' not found"