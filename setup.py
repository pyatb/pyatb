import os
from pathlib import Path
from glob import glob
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import pybind11
from pybind11.setup_helpers import Pybind11Extension

compiler = 'icpc'
mkl_include_dir = None
mkl_library_dir = None
eigen_include_dir = None

include_dirs = []
library_dirs = []
libraries = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'iomp5', 'm', 'dl', 'pthread']
runtime_library_dirs = []

temp_path = Path('siteconfig.py').expanduser()
exec(temp_path.read_text())

include_dirs.append(mkl_include_dir)
include_dirs.append(eigen_include_dir)
include_dirs += [
    os.path.join("src", "cpp", "core"),
    os.path.join("src", "cpp", "interface_python"),
]

library_dirs.append(mkl_library_dir)
runtime_library_dirs.append(mkl_library_dir)

extra_compile_args = []
extra_link_args = []

if compiler == 'g++':
    extra_compile_args += ['-fopenmp']
elif compiler == 'icpc':
    extra_compile_args += ['-qopenmp']

sources = sorted(glob("src/cpp/core/*.cpp")) + sorted(glob("src/cpp/interface_python/*.cpp"))

define_macros = []
undef_macros = []

os.environ["CXX"] = compiler
os.environ["CC"] = compiler
os.environ['CFLAGS'] = "-O2 -std=c++11 -fPIC -Wall -shared"

ext_modules = [
    Pybind11Extension(
        name="pyatb.interface_python",
        sources=sources,
        include_dirs=include_dirs,
        define_macros=define_macros,
        undef_macros=undef_macros,
        library_dirs=library_dirs,
        libraries=libraries,
        runtime_library_dirs=runtime_library_dirs,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language = "c++",
    ),
]

setup(
    name='pyatb',
    version="1.0.0",
    cmdclass={'build_ext': build_ext},
    setup_requires=['numpy', 'pybind11', 'setuptools>=18.0'],
    license='GPL v3.0',
    description='This is the pyatb module.',
    long_description='None',
    author='PYATB Developer',
    author_email='jingan@mail.ustc.edu.com',
    url='https://github.com/pyatb/pyatb',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=['numpy', 'scipy', 'matplotlib'],
    ext_modules=ext_modules,
    entry_points={'console_scripts': ['pyatb = pyatb.main:main']}
)
