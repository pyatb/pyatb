#!/bin/bash

# 加载ABACUS环境
source ~/softwares/abacus-develop-LTSv3.10.0/toolchain/abacus_env.sh

# 运行abacus
export OMP_NUM_THREADS=1 # 设置OpenMP线程数为1
num_mpi=8 # 设置MPI进程数
mpirun -n $num_mpi abacus > job.log 2> job.err