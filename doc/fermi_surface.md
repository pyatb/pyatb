# Fermi surface

## introduction

Fermi surface is a module calculating the iso-energy surface of a given energy. Normally, if the given energy should be Fermi energy, then it would plot Fermi surface.


## example

An example (refer to folder `example/Cu`) of calculating the spectral function of the Cu is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin                  1
    package                ABACUS
    fermi_energy         14.83967742042727
    #fermi_energy           Auto
    #fermi_energy          7
    fermi_energy_unit      eV
    HR_route               data-HR-sparse_SPIN0.csr
    SR_route               data-SR-sparse_SPIN0.csr
    HR_unit                Ry
}

LATTICE
{
    lattice_constant       6.91640 
    lattice_constant_unit  Bohr
    lattice_vector
    0.50	0.50	0.00
    0.50	0.00	0.50
    0.00	0.50	0.50
}
FERMI_SURFACE
{ 
        bar            0.01
        kpoint_mode     mp
        k_start        0 0 0
        k_vect1        1 0.0 0.0
        k_vect2        0.0 1 0.0
        k_vect3        0.0 0.0 1
        mp_grid        50 50 50
}

```
`bar`: The max tolerable error bar for the Fermi surface

`energy`: The given energy. The default value would be the Fermi energy of this system.

For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be two files in the `Out/Fermi_Surface` folder, namely `fermi_surface_kpt.dat` and `plot_fermi_surface.py`,  corresponding to the k-point found on the Fermi surface and a plotting script.
