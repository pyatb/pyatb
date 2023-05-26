# Installation

## Download

The latest version of the PYATB package can be obtained from the GitHub repository:

```shell
git clone https://github.com/pyatb/pyatb.git
```

## Prerequisites

Currently, PYATB is only supported on Linux systems and has not been tested on Windows or Mac systems. To use PYATB, the following prerequisites must be installed:

- Python 3.7 or newer
- C++ compiler
- Intel MKL 
- [pybind11](https://pybind11.readthedocs.io/en/stable/index.html)
- [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [mpi4py](https://mpi4py.readthedocs.io/en/stable/mpi4py.html)
- [Matplotlib](https://matplotlib.org/)


## Install

1. Before installing PYATB, you need to install the [pybind11](https://pybind11.readthedocs.io/en/stable/index.html) module and [mpi4py](https://mpi4py.readthedocs.io/en/stable/mpi4py.html) module:

    ```shell
    pip install pybind11
    ```

    and

    ```shell
    pip install mpi4py
    ```


2. You can set the addresses of the C++ compiler (non-MPI version), MKL library, and Eigen library through the `siteconfig.py` file, and modify the relevant variables to point to the correct paths for your system. 

   - To set the C++ compiler, you can modify the variable "compiler" in the `siteconfig.py`. For example, you can set `compiler = 'icpc'`.

   - To set the MKL library, you can modify the variable "mkl_library_dir" and "mkl_include_dir" in the `siteconfig.py`.

   - To set the Eigen library, you can modify the variable "eigen_include_dir" in the `siteconfig.py`. Eigen is a C++ template library for linear algebra. You can find and download it on the official website [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).

3. After completing the above preparatory work, you can install PYATB with the following command:

    ```shell
    pip install ./
    ```

4. After completing the installation process, you can access the `pyatb` executable and corresponding module, which can be imported using the `import pyatb` command.

## Related to the Intel MKL library

If you encounter the following problem when running pyatb using the intel oneapi MKL library under Anaconda:

```shell
undefined symbol: mkl_sparse_optimize_bsr_trsm_i8
```

This problem may be caused by a conflict between the mkl library in Anaconda and the mkl library in intel oneapi. You can solve it by using the following method:

```shell
export LD_PRELOAD=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:\
/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx2.so.2:\
/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:\
/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:\
/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:\
/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so
```

Note that you need to replace the MKL path with your own environment. For `libiomp5.so`, its location may sometimes be in the directory where the intel compiler is located.

Another solution is to use the MKL library in the Anaconda virtual environment. If you cannot find it, you can install the MKL library and its include files using the following commands:

```shell
conda install -c conda-forge mkl
conda install -c conda-forge mkl-include
```