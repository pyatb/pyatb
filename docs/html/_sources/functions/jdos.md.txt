# JDOS

## introduction

Joint density of states (JDOS) is used to describe the density of states of electrons excited from the valence band to the conduction band, which relates to the absorption spectrum and the dielectric function of system.

The implementation of JDOS per crystal cell is given by

$$
D_{joint}(\omega) = \frac{V_{c}}{\hbar} \int \frac{d^3 k}{(2\pi)^3} \sum_{n,m} f_{nm}\delta(\omega_{mn}-\omega).
$$

where $V_c$ is the cell volume, $f_{nm} = f_n - f_m$ and $\hbar\omega_{mn} = E_{m} - E_{n}$ are differences between occupation factors and band energies, respectively.

Currently JDOS is only used to calculate insulators and semiconductors.

## example

An example of calculating the JDOS of a perovskite CsPbI$_3$ is presented here (refer to folder `examples/CsPbI3`). 

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        5.508975304340945
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            data-rR-sparse.csr
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant        1.8897261258369282
    lattice_constant_unit   Bohr
    lattice_vector
	6.2894000000      0.0000000000      0.0000000000
	0.0000000000      6.2894000000      0.0000000000
	0.0000000000      0.0000000000      6.2894000000
}

JDOS
{
    occ_band      37
    omega         0.5  10
    domega        0.01
    eta           0.2
    grid          30 30 30
}
```

`occ_band`: Specifies the occupied energy band of the system. Currently, only insulator or semiconductor materials can be calculated.

`omega`: Specifies the photon energy, the unit is eV.

`domega`: Specifies the energy interval of the omega.

`eta`: Specify the parameters of Gaussian smearing.

`grid`: Specifies the uniform k-point grid used to calculate the JDOS.

Once the task calculation is finished, two files will be generated in the `Out/JDOS` folder: `JDOS.dat` and `plot_jdos.py`. The former contains the JDOS data, while the latter is a script used for plotting the JDOS.
