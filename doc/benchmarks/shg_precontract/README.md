# SHG three-band precontraction (method=0)

This changes the order of the three-band summation, not the formula. It preserves
all intermediate bands, the existing branch-specific energy cutoffs, occupations,
complex resonance factors, and tensor ordering. The mixed/two-band expressions,
method=1, Python interface, spin factors, and units are unchanged.

## Algebra

For one tensor component, let

```
T(n,m,p) = r_a(n,m) [r_b(m,p) r_c(p,n) + r_c(m,p) r_b(p,n)]
d(n,m,p) = 2 E_p - E_n - E_m
z = omega + i eta
```

The original three branches contribute

```
T/d * [ M2 (f_n-f_m)/(E_m-E_n-2z)
      + Mpn (f_p-f_n)/(2(E_p-E_n-z))
      + Mmp (f_m-f_p)/(2(E_m-E_p-z)) ].
```

Here M2, Mpn, Mmp are the original, distinct cutoff masks. Accumulate the first
coefficient into C2[n,m], the second into C1[n,p], and the third into C1[p,m].
After the triple sum, apply `C1[i,j]/(E_j-E_i-z) + C2[i,j]/(E_j-E_i-2z)`.
Only pairs written by at least one active branch enter the frequency scan.

The regular three-band part changes from O(Nomega Nb^3) to
O(Nb^3 + Nomega Nb^2), with two complex 27-component arrays of size Nb^2
(864 Nb^2 bytes on a platform with 16-byte complex doubles, plus the pair mask).
This does not reduce the cost of interpolation or the existing mixed terms.

## Ill-conditioned static denominators

Triples with |d| < 1e-2 eV retain the original direct frequency-dependent
expressions. If any active triple has |d| < 1e-7 eV, the entire k point retains
the original summation order. These conservative algorithm-selection thresholds
introduce no new broadening or denominator regularization. Worst-case scaling
therefore remains that of the original kernel. Summation-order differences can
still be amplified by cancellation; the checks below are numerical evidence,
not a universal error bound.

The original exact-zero static denominator can produce nonfinite values. This
patch preserves that behavior; resolving it would be a separate physics change.

## Standalone regression check

With Eigen headers available (for example after initializing the Eigen submodule):

```bash
python tests/check_shg_precontract.py --eigen /path/to/eigen --cxx g++
```

The check compiles the actual current arithmetic body and a frozen direct-sum
reference from commit `80f7c2d89aaa8a736e4145f67a6e1d30c783be2d`. Twelve cases
cover both sides of the 1e-2 eV threshold, a near-zero midpoint, an exact singular
midpoint, same-occupation degeneracies, cutoff edges, generic matrices, and
larger band/frequency counts. Eleven finite cases pass a max-norm relative error
bound of 1e-10. One case checks that both versions remain nonfinite. This is not
an assertion that a nonfinite result is physically valid. No timing assertion is
part of the test.

## Measured performance

Baseline: the commit above. Both builds used GCC 9.3.1, Eigen 3.4.0, identical
`-O3 -std=c++17 -fopenmp` flags, and no fast-math. Hardware: Intel Xeon Platinum
8481C, SMT disabled, 28 distinct physical cores on one socket of an exclusively
allocated node. Python 3.10.20, NumPy 2.2.6, SciPy 1.15.3; MPI runs used
OpenMPI 5.0.10 and mpi4py 4.1.2.

A user-supplied 38-band GaAs ABACUS model was used, with a shifted 20^3 mesh
(8000 points), 0.05–3 eV photon energies, eta=0.05 eV, and mu=7.6331377499 eV.
These inputs are not represented as an official PYATB reference dataset.
Four samples per version were collected in ABBA order twice; the table reports
medians of the timed raw C++ response API. MPI uses strided k partitions and
includes SUM reduction and a final barrier in its response timing.

| Parallel mode | Frequencies | Original (s) | Precontracted (s) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| OpenMP, 28 threads | 600 | 55.734 | 9.449 | 5.90x |
| MPI, 28 ranks x 1 thread | 600 | 43.169 | 9.787 | 4.41x |
| OpenMP, 28 threads | 1200 | 107.761 | 14.284 | 7.54x |
| MPI, 28 ranks x 1 thread | 1200 | 81.867 | 13.768 | 5.95x |

All 16 paired runs passed full 27-component complex-array comparisons at
rtol=1e-10, atol=1e-12; the largest component-peak-normalized difference was
1.42e-14. Official and optimized outputs were also compared across parallel
modes for both frequency counts. This is equivalence within the original
finite-band numerical model, not proof of physical convergence.

`results.csv` records timing ranges and workflow timings as well.
`timing_samples.json` contains all paired response times and historical error
metrics (whose original gate used a scale-dependent absolute tolerance).
`validation.json` records the subsequent strict atol=1e-12 checks cited above. OpenMP
workflow timing includes imports, model loading, warmup and export inside the
Python process; MPI workflow timing includes the launcher lifetime. Their
workflow columns therefore have different boundaries and should only be
compared within a mode. The headline table excludes warmup/imports/model IO.
The prior timing campaign used site-specific launch/binding drivers; the
portable runner below reproduces the response calculation and timing boundary,
not those workflow timing columns. No raw material inputs are bundled here.

## Reproducing with supplied GaAs inputs

Input SHA256 (files under `--inputs`):

```
6d3bd27650d0e0ed5921fcebd03c872b5006307d4a1cc3e78adafd7461c06fef  data-HR-sparse_SPIN0.csr
8f12117ec899c919d749fc7ec96c193e0e59948306758a1967a90d19a39b4a81  data-SR-sparse_SPIN0.csr
efe3c3892c0dcdc3b8f356bdd6b66736946ca9e009c1002b6efc97c85c8f594d  data-rR-sparse.csr
```

Build the baseline and candidate in separate checkouts with identical compiler,
Eigen and dependency versions. Point `--source` to each checkout's `src` folder
containing its built extension. Use fresh processes in baseline, candidate,
candidate, baseline order, repeated twice; do not reuse output directories.
Run inside an allocation with 28 physical cores. Restore/verify the full
allocation affinity before every process launch; OpenMP runtime initialization
can pin a launcher thread and accidentally restrict its children to one CPU.

```bash
export OMP_NUM_THREADS=28 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
python doc/benchmarks/shg_precontract/benchmark.py \
  --source /path/to/baseline/src --inputs /path/to/GaAs \
  --grid 20 --omega-num 1200 --output /path/to/new/baseline-run1
```

For MPI, set OMP_NUM_THREADS=1, OMP_PROC_BIND=false, BLAS threads=1; use the
cluster-supported MPI launcher with 28 ranks bound to 28 distinct physical
cores and append `--mpi`. Repeat for both checkouts and for 600 frequencies.
MPI is a benchmark harness partition of the unchanged C++ API; this patch does
not add a new MPI implementation to the solver.

Compare the saved `chi_raw` arrays with `numpy.testing.assert_allclose` using
rtol=1e-10, atol=1e-12 and require identical `omega_eV`/`kpoints` arrays. Do not
apply additional wrapper prefactors before this arithmetic equivalence check.
