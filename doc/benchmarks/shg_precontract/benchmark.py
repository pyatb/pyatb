"""Time the raw C++ method=0 API for a supplied 38-band GaAs ABACUS model.

Run each checkout in a fresh process. This excludes imports/model loading/warmup
from response_seconds. Under MPI it includes reduction and uses max-rank time.
Set and verify CPU affinity externally, before importing NumPy or launching MPI.
The model is user-supplied; its reference hashes are in README.md.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='checkout/src')
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='new directory')
    parser.add_argument('--grid', type=int, default=20)
    parser.add_argument('--omega-num', type=int, default=1200)
    parser.add_argument('--mpi', action='store_true')
    args = parser.parse_args()
    if args.grid < 1 or args.omega_num < 2:
        parser.error('require grid >= 1 and omega-num >= 2')
    if args.output.exists():
        parser.error('--output must not exist')
    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None
    sys.path.insert(0, str(args.source.resolve()))
    os.environ['PYATB_SKIP_RUNTIME_INIT'] = '1'
    import numpy as np
    from pyatb.init_tb import init_tb
    import pyatb.interface_python as extension
    if not Path(extension.__file__).resolve().is_relative_to(args.source.resolve()):
        raise RuntimeError('Imported extension does not belong to --source')
    if args.mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank, ranks = comm.Get_rank(), comm.Get_size()
    else:
        comm, rank, ranks = None, 0, 1
    model = init_tb(
        package='ABACUS', nspin=1, lattice_constant=1.0,
        lattice_vector=np.array([[0, 2.765, 2.765], [2.765, 0, 2.765], [2.765, 2.765, 0]]),
        HR_route=str(args.inputs / 'data-HR-sparse_SPIN0.csr'), HR_unit='Ry',
        SR_route=str(args.inputs / 'data-SR-sparse_SPIN0.csr'), need_rR=True,
        rR_route=str(args.inputs / 'data-rR-sparse.csr'), rR_unit='Bohr', isSparse=False,
    )
    grid = np.indices((args.grid,) * 3).reshape(3, -1).T
    points = np.ascontiguousarray((grid + np.array([.17, .23, .31])) / args.grid)
    local = np.ascontiguousarray(points[rank::ranks])
    if len(local) == 0:
        raise ValueError('Use at most one MPI rank per k point')
    omega = np.linspace(.05, 3., args.omega_num)
    model.tb_solver.get_second_harmonic(0, .05, 3, .2, .2, 7.6331377499, 1, local[:1])
    if comm is not None:
        comm.Barrier()
    begin = time.perf_counter()
    part = model.tb_solver.get_second_harmonic(
        0, .05, len(omega), float(omega[1] - omega[0]), float(omega[0]),
        7.6331377499, len(points), local,
    )
    if comm is not None:
        result = np.empty_like(part) if rank == 0 else None
        comm.Reduce(part, result, op=MPI.SUM, root=0)
        comm.Barrier()
        elapsed = comm.allreduce(time.perf_counter() - begin, op=MPI.MAX)
        bindings = comm.gather(affinity, root=0)
    else:
        result, elapsed, bindings = part, time.perf_counter() - begin, [affinity]
    if rank == 0:
        if not np.isfinite(result).all():
            raise RuntimeError('Nonfinite response')
        args.output.mkdir(parents=True, exist_ok=False)
        np.savez_compressed(args.output / 'tensor.npz', chi_raw=result, omega_eV=omega, kpoints=points)
        metrics = dict(
            nk=len(points), nomega=len(omega), nbands=int(model.basis_num),
            ranks=ranks, response_seconds=elapsed, allocation_affinities=bindings,
            extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
            normalization='raw C++ sum; no Python prefactor or spin multiplication',
        )
        (args.output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
        print(json.dumps(metrics))


if __name__ == '__main__':
    main()
