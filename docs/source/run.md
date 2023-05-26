# Run

PYATB supports mixed parallelism of OpenMP and MPI, and you need to determine the number of threads and processes to run depending on the actual configuration of your computer. 

For example, set the number of threads to 2,
```shell
export OMP_NUM_THREADS=2
```
and then use 6 processes to run PYATB,
```shell
mpirun -n 6 pyatb
```