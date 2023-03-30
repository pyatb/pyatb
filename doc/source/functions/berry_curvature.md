# BERRY_CURVATURE

## introduction

Berry curvature is of fundamental importance for understanding some basic properties of solid materials and is essential for the description of the dynamics of Bloch electrons.

The Berry curvature of a single energy band is defined as follows:

$$
\Omega_n(\mathbf{k}) = \nabla \times \mathbf{A}_{n}(\mathbf{k}),
$$ 

where Berry phase $\mathbf{A}_{n}(\mathbf{k}) = i \langle u_{n\mathbf{k}}|\nabla_{\mathbf{k}}|u_{n\mathbf{k}}\rangle$, $|u_{n\mathbf{}k}\rangle$ is the periodic part of the Bloch wave function.

We calculated the Berry curvature:

$$
\Omega_{\alpha\beta}(\mathbf{k}) = \sum_{n} f_n(\mathbf{k}) \Omega_{n, \alpha\beta}(\mathbf{k}),
$$ 

where $f_n$ is the Fermi occupation function.

Detailed descriptions can be found in [Calculation of Berry curvature using non-orthogonal atomic orbitals](https://doi.org/10.1088/1361-648X/ac05e5).

## example

An example (refer to folder `example/Fe`) of calculating the Berry curvature of the fcc-Fe is given here.

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin               4
    package             ABACUS
    fermi_energy        18.18839115931923
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            data-rR-sparse.csr
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant        5.4235
    lattice_constant_unit   Bohr
    lattice_vector
     0.5  0.5  0.5
    -0.5  0.5  0.5
    -0.5 -0.5  0.5
}

BERRY_CURVATURE
{
    method                  0
    kpoint_mode             line
    kpoint_num              10
    high_symmetry_kpoint
    0.0   0.0    0.0   100 # G
    0.5  -0.5   -0.5   100 # H
    0.75  0.25  -0.25  100 # P
    0.5   0.0   -0.5   100 # N
    0.0   0.0    0.0   100 # G
    0.5   0.5    0.5   100 # H
    0.5   0.0    0.0   100 # N
    0.0   0.0    0.0   100 # G
    0.75  0.25  -0.25  100 # P
    0.5   0.0    0.0   1   # N
}
```

`method`: Method for calculating Berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.

`occ_band`: The number of occupied energy bands of an insulator. When this value is not set, it will be determined according to the Fermi energy.

For the k point setting of this function, please refer to the `kpoint_mode` module.

Upon completion of the task calculation, two files will be generated in the `Out/Berry_Curvature` directory: `kpt.dat` and `berry_curvature.dat`. These files respectively contain the k-point coordinates and the total Berry curvature for each k-point.
