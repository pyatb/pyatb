# COHP

The COHP function evaluates atom-pair crystal orbital Hamilton population
spectra from the ABACUS Hamiltonian matrix, overlap matrix and NAO
eigenvectors already used by PyATB.

The implemented value for one band and one k point is

```text
sum_ij Re[ C_i,n(k)^* M_ij(k) C_j,n(k) ]
```

where `M` is `H(k)` for COHP and `S(k)` for COOP. The selected orbital sets
belong to two atom indices in the ABACUS `STRU` file. The spectrum is broadened
on an energy grid with the same Gaussian smearing convention used by the PyATB
PDOS/JDOS modules.

## Input block

```text
COHP
{
    stru_file          STRU
    atom_i_index      1
    atom_j_index      2
    atom_i_orbs       all
    atom_j_orbs       2p
    method            COHP
    spin              sum
    e_range           -8.0 8.0
    de                0.02
    sigma             0.08
    invert            1
    shift_to_efermi   1
    output_prefix     COHP
    input_file        Input
    orbital_dir       basis
    kpoint_mode       mp
    mp_grid           8 8 8
}
```

`atom_i_orbs` and `atom_j_orbs` accept `all`, shell selectors such as `s`,
`p`, `d`, shell-like aliases such as `2s`, `2p`, `3d`, or comma/space
separated global orbital indices. The numeric prefix is currently treated only
as an angular-momentum alias: for example, `2p` is equivalent to `p`. It does
not select a distinct radial shell. Atom indices are 1-based and follow the
atom order in `STRU`.

`orbital_dir` points to the directory containing ABACUS numerical orbital
files. If `orbital_dir` is omitted, COHP tries to read `orbital_dir` from
`input_file`; when the COHP block omits `input_file`, PyATB passes the running
`Input` file path. If neither source gives an orbital directory, COHP falls
back to the directory containing `STRU`.

## Output

The output directory is `Out/COHP`. The main spectrum is written to
`COHP.dat` by default, with columns `energy_eV` and the signed method label (default `-COHP`). If
`shift_to_efermi` is set to `1`, the energy axis is `E - E_F`.

`COHP.meta.json` records the selected atoms, selected global orbital indices,
orbital shell mapping and output file names. `plot_cohp.py` is generated in the
same directory and can be run explicitly to create `cohp.pdf`.

## Notes

COHP currently supports `nspin = 1` and collinear `nspin = 2`. For spin-polarized
data, set `spin` to `sum`, `up`, or `down`.

## Pair, sign and integration conventions

The two resolved **global orbital sets**, rather than the atom IDs, determine
pair normalization: identical sets use factor 1, disjoint sets use factor 2
(the sum of both Hermitian directions). This also applies to disjoint channels
on the same atom. Partially overlapping or empty selections raise `ValueError`;
split overlapping queries into identical/disjoint blocks explicitly.

`invert=1` (the default) writes `-COHP` or `-COOP`; `invert=0` writes the raw
population. The data header, metadata `spectrum_label`, and generated plot agree.
Existing numeric columns and filenames are retained. COHP spectra have units
eV/eV (dimensionless), COOP spectra 1/eV.

`COHP.meta.json` additionally contains `occupied_integral`:

```json
{
  "method": "ICOHP",
  "occupation": "E <= EF (zero temperature)",
  "unit": "eV",
  "raw_value": -1.0,
  "value": 1.0,
  "label": "-ICOHP",
  "spin_degeneracy": 2.0
}
```

The numbers above are illustrative. This is the discrete sum over eigenstates
with **absolute energy E <= fermi_energy**, including normalized k weights and
the selected spin channel(s), before any broadening. `raw_value` preserves the
Hamilton/overlap population sign; `value` follows `invert`. ICOOP uses the same
object with method `ICOOP` and unit `dimensionless`. These are zero-temperature
cutoff sums, not finite-temperature occupation-weighted bond energies; SCF
fractional occupations are not an input. They are independent of `e_range`,
`de`, `sigma`, and shifting the displayed energy axis. A line-path k sample
still yields its weighted sum, not a Brillouin-zone integral.

K weights default to uniform and explicit non-negative finite weights are
normalized to sum to one. An explicit `kpoint_weights` keyword without values
is an error; omit the keyword to request uniform weights. `nspin=1, spin=sum` has degeneracy 2. For `nspin=2`,
`up`/`down` select one channel and `sum` adds both, with no extra factor 2.

The Gaussian retains PyATB's convention
`exp(-(dE/sigma)^2)/(sqrt(pi)*sigma)` (standard deviation `sigma/sqrt(2)`).
Contributions within 10 sigma are sampled even when the state is at or outside
the plotting window. The finite plotted portion is **not renormalized** and its
numerical integral need not equal the discrete occupied integral. Grid samples
start at `e_range[0]` with spacing `de`, ending at the last sample at or below
`e_range[1]` (up to floating-point rounding). Metadata records requested and
sampled ranges and the broadening convention. `fermi_energy` and range endpoints
must be finite; `de` and `sigma` must be finite and positive, and the range must
span at least one grid step.
