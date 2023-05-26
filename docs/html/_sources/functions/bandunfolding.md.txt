# BANDUNFOLDING

## introduction

Band unfolding is a method of projecting the wave function of the supercell to the coupled $\mathbf{k}$ points in the original unit cell to obtain the spectral function.

The relationships between the lattice vectors of the large cell (A) and  primitive cell (a) are given by

$$
\begin{bmatrix}
    A_1 \\
    A_1 \\
    A_3
\end{bmatrix} = 
\begin{bmatrix}
    m_{11} & m_{12} & m_{13} \\
    m_{21} & m_{22} & m_{23} \\
    m_{31} & m_{32} & m_{33} \\
\end{bmatrix}
\begin{bmatrix}
    a_1 \\
    a_2 \\
    a_3
\end{bmatrix}
$$

Detailed descriptions can be found in [First-principles calculations of the surface states of doped and alloyed topological materials via band unfolding method](https://doi.org/10.1016/j.commatsci.2022.111656).

## example

Here is an example located in the `examples/NV` folder that showcases how to calculate the spectral function of the nitrogen-vacancy (NV) center in diamond.

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                           1   
    package                         ABACUS
    fermi_energy                    15.46412271260700
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
    7.13366 0 0 #latvec1
    0 7.13366 0 #latvec2
    0 0 7.13366 #latvec3
}

BANDUNFOLDING
{
    stru_file                       STRU
    ecut                            60
    band_range                      10 250
    m_matrix                        -2 2 2 2 -2 2 2 2 -2
    kpoint_mode                     line
    kpoint_num                      5
    high_symmetry_kpoint
    0.500000  0.000000  0.500000 300  # X
    0.500000  0.250000  0.750000 300  # W
    0.500000  0.500000  0.500000 300  # L
    0.000000  0.000000  0.000000 300  # Gamma
    0.500000  0.000000  0.500000 1    # X
}
```

`stru_file`: The structure file name of the supercell. This file indicates the crystal structure of the supercell and the corresponding orbital file. Make sure that both the structure file and the orbital file exist.

`ecut`: Determine the number of projections to the plane wave basis group, the energy unit is Ry.

`band_range`: Specify the energy band range of band unfolding.

`m_matrix`: Transformation matrix of supercell and primitive cell lattice vector

The k-point setting is determined according to the structure of the original cell, not the k-point in the supercell. For the k point setting of this function, please refer to the `kpoint_mode` module.

Once the task calculation is finished, three files will be generated in the `Out/Bandunfolding` folder. These files include `kpt.dat` and `spectral_weight.dat`, representing the k-point and spectral function of the original cell, respectively. Additionally, a drawing script called `plot_unfold.py` will also be included.

**NOTE**: When calculating a slab material, it is necessary to remove the thickness of the vacuum layer and make necessary modifications to the crystal vector and atomic positions specified in the `LATTICE` and structure files. Failing to make these modifications accurately can result in inaccurate results.
