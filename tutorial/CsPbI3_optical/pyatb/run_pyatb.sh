#!/bin/bash

# 加载conda环境
source ~/software/miniconda3/bin/activate
conda activate pyatb

# 运行pyatb
export OMP_NUM_THREADS=1 # 设置OpenMP线程数为1
num_mpi=8 # 设置MPI进程数
mpirun -n $num_mpi pyatb > job.log 2> job.err