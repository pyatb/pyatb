#!/bin/bash

# 构建自定义 Docker 镜像
docker build -t pyatb-custom-manylinux .

# 使用自定义镜像进行构建
CIBW_MANYLINUX_X86_64_IMAGE=pyatb-custom-manylinux CIBW_BUILD="cp313-manylinux_x86_64" python -m cibuildwheel --output-dir dist
