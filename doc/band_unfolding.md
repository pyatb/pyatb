# Band unfolding

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


## example

An example (refer to folder `example/MnBi2Te4-afm`) of calculating the spectral function of the AFM MnBi$_2$Te$_4$ slab is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin               4
    package             ABACUS
    fermi_energy        9.540272009417667
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    HR_unit             Ry
    rR_unit             Bohr
    max_kpoint_num      8000
}

LATTICE
{
    lattice_constant        1.8897162
    lattice_constant_unit   Bohr
    lattice_vector
    4.3337001801 0 0 
    2.16685009 3.7530944483 0 
    4.3339143734 2.50218663 27.2734777344 
}

BANDUNFOLDING
{
    stru_file               STRU
    ecut                    40
    band_range              176 276
    m_matrix                1 0 0 0 1 0 0 0 2
    kpoint_mode             line
    kpoint_num              5
    high_symmetry_kpoint
    0.0 0.0 0.0 100 #Gamma
    0.5 0.0 0.5 100 #Z
    0.5 0.5 0.0 100 #F
    0.0 0.0 0.0 100 #Gamma 
    0.5 0.0 0.0 1   #L
}
```
`stru_file`: The structure file name of the supercell. This file indicates the crystal structure of the supercell and the corresponding orbital file. Make sure that both the structure file and the orbital file exist.

`ecut`: Determine the number of projections to the plane wave basis group, the energy unit is eV.

`band_range`: Specify the energy band range of band unfolding.

`m_matrix`: Transformation matrix of supercell and primitive cell lattice vector

The k-point setting is determined according to the structure of the original cell, not the k-point in the supercell. For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be three files in the `Out/Bandunfolding` folder, namely `kpt.dat` and `spectral_weight.dat`, `plot_unfold.py`, corresponding to the k-point and spectral function of the original cell, the drawing script.
