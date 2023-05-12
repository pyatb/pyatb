import os
# from ctypes.util import find_library
import sysconfig
from glob import glob
import pybind11
from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup

include_dirs = [
    os.path.join("src", "cpp", "core"),
    os.path.join("src", "cpp", "interface_python"),
    sysconfig.get_path("include"),
    pybind11.get_include(),
    os.path.join(sysconfig.get_path("purelib"), "cmeel.prefix", "include", "eigen3")
]

sources = sorted(glob("src/cpp/core/*.cpp")) + sorted(glob("src/cpp/interface_python/*.cpp"))

define_macros = []
undef_macros = []

os.environ["LDFLAGS"] = sysconfig.get_config_var("LDFLAGS")

libraries = []
# libraries.append(find_library('libmkl_intel_lp64.so'))
# libraries.append(find_library('libmkl_gnu_thread.so'))
# libraries.append(find_library('libmkl_core.so'))

extra_compile_args = ["-fopenmp", "-lmkl_intel_lp64", "-lmkl_gnu_thread", "-lmkl_core", "-lgomp", "-lm", "-lpthread",]

ext_modules = [
    Pybind11Extension(
        name="pyatb.interface_python",
        sources=sources,
        include_dirs=include_dirs,
        define_macros=define_macros,
        undef_macros=undef_macros,
        # library_dirs=library_dirs,
        libraries=libraries,
        # runtime_library_dirs=runtime_library_dirs,
        extra_compile_args=extra_compile_args,
        # extra_link_args=extra_link_args,
        language = "c++",
    ),
]

setup(
    ext_modules=ext_modules,
)
