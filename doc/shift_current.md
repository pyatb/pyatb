# Shift current

## introduction

The shift current is one of the intrinsic contributions of the bulk photovoltaic effect (BPVE). The photovoltage generated through the shift current mechanism can exceed the energy gap of the system and is not constrained by the Shockley-Queisser limit.

$$
J^{a} = 2 \sigma^{abc}(0; \omega, -\omega) E_{b}(\omega) E_{c}(-\omega)
$$

Where $\sigma^{abc}(0; \omega, -\omega)$ is the shift current conductivity.

## example

An example (refer to folder `example/GeS`) of calculating the shift current of the GeS monolayer is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        0.4641700721815871
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            new-data-rR-tr_SPIN1
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant        1.0
    lattice_constant_unit   Angstrom
    lattice_vector
    15.00 0.00 0.00
    0.00 3.66 0.00
    0.00 0.00 4.47
}


SHIFT_CURRENT
{
    occ_band      20
    omega         0   4
    domega        0.01
    eta           0.01
    grid          1 1000 1000
    method        1
}
```

`occ_band`: Specifies the occupied energy band of the system. Currently, only insulator or semiconductor materials can be calculated.

`omega`: Specifies the photon energy, the unit is eV.

`domega`: Specifies the energy interval of the omega.

`eta`: Specify the parameters of Gaussian smearing.

`grid`: Specifies the uniform k-point grid used to calculate the shift current.

`method`: There are currently two different methods to solve shift current.

After completing the task, two main files are generated in the `Out/Shift_current` folder, namely `shift_current.dat` and `plot_shift_current.py`.