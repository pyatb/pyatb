import os
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

# 用户修改的参数
CXX = "icpc"
LAPACK_DIR = "/opt/intel/oneapi/mkl/2022.0.2"


# 一般不需要修改
os.environ["CXX"] = CXX
os.environ["CC"] = CXX
os.environ['CFLAGS'] = "-O2 -std=c++11 -fPIC -Wall -shared"
LAPACK_INCLUDE_DIR = LAPACK_DIR + "/include"
LAPACK_LIB_DIR = LAPACK_DIR + "/lib/intel64"
LAPACK_LIB = "-L" + LAPACK_LIB_DIR + " -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -Wl,-rpath=" + LAPACK_LIB_DIR
MY_INCLUDE = "cpp/include"
sources_dir = "cpp/src/core"
sources = [
    "band_structure_solver.cpp",
    "bandunfolding_solver.cpp",
    "base_data.cpp",
    "berry_connection_solver.cpp",
    "berry_curvature_solver.cpp",
    "berry_phase_solver.cpp",
    "cell_atom.cpp",
    "linear_response.cpp",
    "math_integral.cpp",
    "math_sphbes.cpp",
    "optical_conductivity_solver.cpp",
    "shift_current_solver.cpp",
    "tools.cpp",
    "velocity_solver.cpp",
    "xr_operation.cpp"
]

for i, s in enumerate(sources):
    sources[i] = os.path.join(sources_dir, s)

include_dirs = [sources_dir, MY_INCLUDE]
extra_compile_args = []
if CXX == "g++":
    extra_compile_args += ["-fopenmp"]
else:
    extra_compile_args += ["-qopenmp"]

extra_link_args = extra_compile_args + [LAPACK_LIB]
define_macros = []
undef_macros = []

extension = Extension('pyatb.interface_python',
                      include_dirs=include_dirs,
                      sources=["cpp/src/interface_python/interface_python.cpp"] + sources,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args,
                      define_macros=define_macros,
                      undef_macros=undef_macros)

setup(name='pyatb',
      version="1.0.0",
      cmdclass={'build_ext': build_ext},
      setup_requires=['numpy', 'setuptools>=18.0'],
      license='GPL v3.0',
      description='This is the pyatb module.',
      long_description='None',
      author='Gan Jin & Hongsheng Pang',
      author_email='jingan@mail.ustc.edu.com',
      url='None',
      packages=find_packages(),
      py_modules=[],
      install_requires=['numpy', 'mpi4py', 'scipy', 'mpi4py'],
      ext_modules=[extension],
      entry_points={'console_scripts': ['pyatb = pyatb.main:main']})



