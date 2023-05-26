# POLARIZATION

## introduction

We calculate the spontaneous polarization of periodic solids by so-called Modern Theory of Polarization, namely Berry phase theory. The electric polarization $\mathbf{P}$ is a modulo a quantum $e\mathbf{R}/V_{c}$ multi-valued function, corresponding to the following implementation equation:

$$
\mathbf{P} = \frac{-e}{(2\pi)^3} \sum_{n}^{occ} \int_{\text{BZ}} \mathbf{A}_{n}(\mathbf{k}) d^3k,
$$

where $\mathbf{A}_{n}(\mathbf{k})$ is the Berry connection of a single band.

## example

Here, we present an example (located in the `examples/PbTiO3` folder) showcasing the calculation of the polarization of PbTiO$_3$.

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        13.38267075814371
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            data-rR-sparse.csr
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant        7.3699
    lattice_constant_unit   Bohr
    lattice_vector
    1.0000000000         0.0000000000         0.0000000000
    0.0000000000         1.0000000000         0.0000000000
    0.0000000000         0.0000000000         1.0000000000
}

POLARIZATION
{
    occ_band      22
    nk1           10
    nk2           10
    nk3           10
    atom_type     3
    stru_file     STRU
    valence_e     14 12 6 
}
```

`occ_band`: The number of occupied energy bands of an insulator.

`nk1`: This refers to the number of samples taken in the $x$ direction of the reciprocal lattice vector $\mathbf{G}$.

`nk2`: This refers to the number of samples taken in the $y$ direction of the reciprocal lattice vector $\mathbf{G}$.

`nk3`: This refers to the number of samples taken in the $z$ direction of the reciprocal lattice vector $\mathbf{G}$.

`stru_file`: Specify the strucutre file. NAOs files are not required.

`atom_type`: The number of element types in the system.

`valence_e`: The number of valence electrons per element.

The parameters `nk1`, `nk2`, and `nk3` correspond to the number of k-points used for the integration in the three lattice directions, and increasing their values leads to a more accurate and convergent result.

After completing the task, `polarization.dat` appears in the `Out/Polarization` folder which contains the electric polarization of the three lattice directions.
