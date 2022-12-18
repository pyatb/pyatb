# AHC

## introduction

The dc anomalous Hall conductivity (AHC) is simply given as the Brillouin zone integral of the Berry curvature of occupying energy bands, 

$$
\sigma_{xy} = -\frac{e^2}{\hbar} \sum_{n}^{occ} \int_{\text{BZ}} \frac{d\mathbf{k}}{(2\pi)^3} f_n(\mathbf{k})\Omega_{n, z}(\mathbf{k}).
$$

## example

An example (refer to folder `example/Fe`) of calculating the AHC of the fcc-Fe is given here.

The `Input` file are:

```text {.line-numbers}
INPUT_PARAMETERS
{
    nspin               4
    package             ABACUS
    fermi_energy        18.18839115931923
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            new-data-rR-tr_SPIN4
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

AHC
{
    integrate_mode          Grid
    integrate_grid          100 100 100
    adaptive_grid           20 20 20
    adaptive_grid_threshold 100  
}
```

`integrate_mode`: Specifies the mode of integration, which can be lattice integration and adaptive integration.

`integrate_grid`: Specifies a uniform grid for grid integration.

`adaptive_grid`: Specifies the grid for adaptive densification.

`adaptive_grid_threshold`: Specifies the cut-off value of adaptive densification, the unit is $\AA^2$.