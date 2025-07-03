import os
import subprocess
from glob import glob
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
try:
    from pybind11.setup_helpers import Pybind11Extension
except:
    from setuptools import Extension as Pybind11Extension

libraries = ['openblas', 'lapacke']

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

# 设置默认的大版本号和是否是 release 版本
ver_default = '1.1.1'
ver_release = True           # 将此设置为 True 以标记为 release 版本

# 获取环境变量中的 commit hash 值
hash_suffix = os.getenv('HASH', '')

# 如果不是 release 版本，且没有从环境变量获取到 hash 值，则尝试从本地 Git 仓库获取
if not ver_release and not hash_suffix:
    try:
        # 获取当前的 Git commit hash
        hash_suffix = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('utf-8').strip()
    except subprocess.CalledProcessError:
        # 如果获取失败，hash_suffix 仍然为空
        hash_suffix = ''

# 拼接版本号
# 如果是 release 版本，使用默认的大版本号
# 如果不是 release 版本，且有 hash_suffix，则生成包含 hash 的版本号
version = f"{ver_default}" if ver_release else f"{ver_default}.dev+{hash_suffix}" if hash_suffix else f"{ver_default}.dev"

# print(f"Building version check HASH: {version}")  # 打印版本号，确保它正确

with open("README.md", encoding="utf-8") as f:
    lines = f.readlines()
    long_description = "".join(lines[1:])

setup(
    name='pyatb',
    version=version,  # 使用动态版本号 with commit hash
    cmdclass={'build_ext': build_ext},
    setup_requires=["setuptools>=42", "wheel", "pybind11"],
    license='GPL v3.0',
    description='This is the pyatb module.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='PYATB Developer',
    author_email='jingan@mail.ustc.edu.cn',
    url='https://github.com/pyatb/pyatb',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    # install_requires=['numpy >= 1.23.5', 'scipy >= 1.9.0', 'mpi4py < 4.0.0', 'matplotlib >= 3.6.0', 'ase'],
    install_requires=['numpy >= 1.17.0', 'scipy >= 1.5.0', 'mpi4py >= 3.1.0', 'matplotlib >= 2.2.2', 'ase'],
    ext_modules=ext_modules,
    entry_points={'console_scripts': [
        'pyatb = pyatb.main:main',
        'pyatb_input = pyatb.easy_use.input_generator:main'
    ]}
)
