from mpi4py import MPI
import numpy as np

# mpi parameters
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

def gather_numpy(a, b, dt):
    return np.r_[a, b]

op_gather_numpy = MPI.Op.Create(gather_numpy, commute=False)

op_sum = MPI.SUM