# ABACUS CSR input compatibility

The ABACUS H/S readers accept both the legacy `Matrix Dimension` / `Matrix number`
layout and the annotated HContainer layout (`number of localized basis`,
`number of Bravais lattice vector R`, `CSR Format`). File names do not select the
parser. The modern layout supports wrapped values/indices/pointers, comments,
empty blocks, and complex `(real,imag)` entries. Truncated arrays, invalid CSR
pointers/indices, and non-finite H/S values raise `ValueError`.

The readers retain the first matrix snapshot, upper-triangular storage, HR unit
conversion (Ry/eV), and SR normalization used by the legacy implementation.
The rR reader remains the legacy-layout reader; the observed develop `rr.csr`
also uses that layout. This does not promise support for an unobserved new rR
layout or every ABACUS output variant.

ABACUS may omit different all-zero lattice-vector blocks in H, S and rR.
The CLI and Python `init_tb` ABACUS paths align all supplied matrices on the
ordered union of their R vectors before setting up the solver. Missing sparse
blocks are zero; S-only, spin-down-only, and rR-only vectors are retained.
The input matrix objects are not modified. Solver interfaces and WANNIER90
loading are unchanged.

## Regression and validation scope

Run `python -m pytest -q tests/test_abacus_csr_reader.py`. The tests cover
legacy/annotated layouts, wrapped/complex data, malformed arrays, non-finite
values, zero blocks, and nspin=1/2/4 handoff and R alignment. Tests use lightweight
MPI/native doubles; they do not require launching ABACUS or a cluster.

The compressed real Si2/LCAO/nspin=1 H/S fixtures were produced with ABACUS
`94576a801` (develop) and `f71921fe8` (LTS). Producer identity and SHA-256 hashes
of the uncompressed inputs are recorded in `tests/fixtures/abacus-csr/provenance.json`.
Both H inputs have 93 R blocks; S has 63 in develop and 93 in LTS, whose extra
blocks are zero. After alignment, the stored H/S values agree exactly.

The real fixtures cover Si/nspin=1. Synthetic nspin=2/4 tests establish input
handling only, not real spin-polarized/SOC calculations or scientific convergence.
