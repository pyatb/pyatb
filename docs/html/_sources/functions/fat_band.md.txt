# FAT_BAND

## introduction

A fat band can provide information about the contributions of specific atomic orbitals or groups of orbitals to the electronic bands of a material at given k points.

## example

In the `examples/Si2` folder, there is an example of how to calculate the fat band of the diamond Si.

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        6.389728305291531
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            data-rR-sparse.csr
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant        1.8897162
    lattice_constant_unit   Bohr
    lattice_vector
    0.000000000000  2.715000000000  2.715000000000
    2.715000000000  0.000000000000  2.715000000000
    2.715000000000  2.715000000000  0.000000000000
}

FAT_BAND
{
    band_range                     1 8
    stru_file                      STRU
    kpoint_mode                    line
    kpoint_num                     5
    high_symmetry_kpoint
    0.50000  0.50000 0.5000 100  # L
    0.00000  0.00000 0.0000 100  # G
    0.50000  0.00000 0.5000 100  # X
    0.37500 -0.37500 0.0000 100  # K
    0.00000  0.00000 0.0000 1    # G
}
```

`band_range`: There are two numbers (separated by spaces) to indicate which bands are selected for projection, counting from 1.

`stru_file`: The structure file name. This file indicates the crystal structure and the corresponding orbital file. Make sure that both the structure file and the orbital file exist.

For the k point setting of this function, please refer to the `kpoint_mode` module.

Once the task calculation is finished, you will find four files in the `Out/Fat_Band` folder. These files include `band.dat` and `pband.dat`, `fatband.xml`, and `plot_fatband.py`. They contain valuable information about the original bands, the coefficients of the bands projected onto each atomic orbital (the number of atomic orbitals is equal to the number of basis sets), an XML formatted file of the projected bands, and a script to visualize the fat band.
