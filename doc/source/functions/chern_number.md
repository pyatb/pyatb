# CHERN_NUMBER

## introduction

Chern number is a topological invariant used to explain the quantized Hall conductivity.  The Chern invariant is the total Berry flux in the 2D Brillouin zone,

$$
n = \frac{1}{2\pi} \oint_{\mathbf{S}} \boldsymbol{\Omega}\cdot d\mathbf{S} 
$$

where $n$ is an integer, $\mathbf{S}$ is any closed 2D manifold, $\boldsymbol{\Omega}$ is total Berry curvature (flux).

To calculate the Chern number, you must first select a closed 2D surface in the Brillouin.

## example

An example (refer to folder `examples/MnBi2Te4-weyl`) of calculating Chern number of the Weyl semimetal MnBi$_2$Te$_4$ is given here.

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                           4
    package                         ABACUS
    fermi_energy                    9.2309138700265265
    fermi_energy_unit               eV
    HR_route                        data-HR-sparse_SPIN0.csr
    SR_route                        data-SR-sparse_SPIN0.csr
    rR_route                        data-rR-sparse.csr
    HR_unit                         Ry
    rR_unit                         Bohr
}

LATTICE
{
    lattice_constant                1.8897162
    lattice_constant_unit           Bohr
    lattice_vector
    4.3773399999000002 0.0000000000000000  0.0000000000000000
    2.1886700000000001 3.7908876409999999  0.0000000000000000
    2.1886700000000001 1.2636292099999999 13.7730333333000008
}

CHERN_NUMBER
{
    method                          0
    occ_band                        109
    integrate_mode                  Grid
    integrate_grid                  100 100 1
    adaptive_grid                   20  20  1
    adaptive_grid_threshold         100
    k_start                         0 0 0
    k_vect1                         1 0 0
    k_vect2                         0 1 0
}
```

`method`: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.

`occ_band`: The number of occupied energy bands of an insulator. When this value is not set, it will be determined according to the Fermi energy.

To determine a plane of k-space requires an origin (`k_start`) and two vectors that are not parallel to each other (`k_vect1`, `k_vect2`).

`k_start`: The origin point coordinates used to describe a Brillouin zone plane.

`k_vect1`: The expansion vector used to describe a Brillouin zone plane.

`k_vect2`: The expansion vector used to describe a Brillouin zone plane.

For k-point integration, please refer to the [Setting of integration] section.

After the task calculation is completed, the `chern_number.dat` file will appear in the `Out/Chern_Num` folder which contains the Chern number specific results.
