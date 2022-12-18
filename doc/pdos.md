# PDOS

# introduction

PDOS is to project the Bloch wave function onto the basis set of atomic orbitals and calculate the contribution of different orbitals to DOS.

# example

An example (refer to folder `example/dos_si2`) of calculating the PDOS of the diamond Si is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        6.585653951974235
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant        10.2
    lattice_constant_unit   Bohr
    lattice_vector
    0.5 0.5 0.0
    0.5 0.0 0.5
    0.0 0.5 0.5
}

PDOS
{
    stru_file     STRU
    e_range       -5.0 17.0
    de            0.01
    sigma         0.07
    kpoint_mode   mp
    mp_grid       12 12 12
}
```

`stru_file`: The structure file name of the supercell. This file indicates the crystal structure of the supercell and the corresponding orbital file. Make sure that both the structure file and the orbital file exist.

`e_range`: Specify the energy range of dos, the unit is eV.

`de`: specifies the energy interval.

`sigma`: Specify the parameters of Gaussian smearing.

The k-point setting is determined according to the structure of the original cell, not the k-point in the supercell. For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be three files in the `Out/PDOS` folder, namely `TDOS.dat` and `PDOS.xml`, `plot_dos.py`. Specify the projected atomic orbital index in the plot script, and then draw the PDOS plot.

