# FERMI_SURFACE

## introduction

The Fermi surface is a function that calculates the iso-energy surface of a given energy level. If the given energy level is the Fermi energy, then it will plot the Fermi surface. 

## example

You can find an example of how to calculate the Fermi surface of Copper in the `examples/Cu` folder. 

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                  1
    package                ABACUS
    fermi_energy           Auto
    fermi_energy_unit      eV
    HR_route               data-HR-sparse_SPIN0.csr
    SR_route               data-SR-sparse_SPIN0.csr
    HR_unit                Ry
    rR_unit                Bohr
    max_kpoint_num         8000
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

FERMI_ENERGY
{
    temperature            0
    electron_num           11
    grid                   50 50 50
    epsilon                1e-4
}

FERMI_SURFACE
{
    bar                    1e-5
    kpoint_mode            mp
    k_start                0 0 0
    k_vect1                1 0 0
    k_vect2                0 1 0
    k_vect3                0 0 1
    mp_grid                50 50 50
}
```

`bar`: The max tolerable error bar for the Fermi surface.

`nbands`: If you know the energy band range where the Fermi energy is located, you can set this parameter to speed up the calculation. There are two numbers in total, indicating the range of the energy band. The default value is `0 0`, that is, all energy bands are considered.

For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be two files in the `Out/Fermi_Surface` folder, namely `fermi_surface_kpt.dat` and `plot_fermi_surface.py`,  corresponding to the k-point found on the Fermi surface and a plotting script.

**NOTE**: To visualize iso-energy surface at specific energy levels, the `fermi_energy` parameter in `INPUT_PARAMETERS` can be customized. 
