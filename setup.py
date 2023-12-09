import os
from glob import glob
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
try:
    from pybind11.setup_helpers import Pybind11Extension
except:
    from setuptools import Extension as Pybind11Extension

libraries = ['lapack', 'blas']

include_dirs = [
    os.path.join("src", "cpp", "core"),
    os.path.join("src", "cpp", "interface_python"),
    os.path.join("eigen"),
]

extra_compile_args = ['-fopenmp']
extra_link_args = ['-lgomp']

sources = sorted(glob("src/cpp/core/*.cpp")) + sorted(glob("src/cpp/interface_python/*.cpp"))

ext_modules = [
    Pybind11Extension(
        name="pyatb.interface_python",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language="c++",
    ),
]

setup(
    name='pyatb',
    version="1.0.0",
    cmdclass={'build_ext': build_ext},
    setup_requires=['pybind11', 'setuptools>=18.0'],
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
