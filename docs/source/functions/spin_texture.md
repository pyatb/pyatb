# SPIN_TEXTURE

## introduction

Spin texture refers to the spatial distribution of electron spins in momentum space. In PYATB, the spin texture is calculated as follows,

$$
\langle \Psi_{n\mathbf{k}} | \hat{\sigma}_i  | \Psi_{n\mathbf{k}} \rangle
= \sum_{\mu,\nu,s,s^{\prime}} C^{*}_{n,\mu s}(\mathbf{k}) S_{\mu\nu, s s^{\prime}}(\mathbf{k}) \hat{\sigma}_{i,s s^{\prime}} C_{n,\nu s^{\prime}}(\mathbf{k}),
$$

where $\hat{\sigma}_i$ are the Pauli matrices, with $i$= $x$, $y$, $z$, and $s$=$\uparrow$, $\downarrow$ is the spin index.

## example

Here, we provide an example of calculating the spin texture of Bi$_2$Se$_3$ (refer to folder `examples/Bi2Se3`). 

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                          4
    package                        ABACUS
    fermi_energy                   9.557219691497478
    fermi_energy_unit              eV
    HR_route                       data-HR-sparse_SPIN0.csr
    SR_route                       data-SR-sparse_SPIN0.csr
    rR_route                       data-rR-sparse.csr
    HR_unit                        Ry
    rR_unit                        Bohr
    max_kpoint_num                 8000
}

LATTICE
{
    lattice_constant               1.8897162
    lattice_constant_unit          Bohr
    lattice_vector
    -2.069  -3.583614  0.000000
     2.069  -3.583614  0.000000
     0.000   2.389075  9.546667
}

SPIN_TEXTURE
{
    nband              78
    kpoint_mode        direct
    kpoint_num         20
    kpoint_direct_coor
    0.010000  0.000000 0.000000
    0.009511  0.003090 0.000000
    0.008090  0.005878 0.000000
    0.005878  0.008090 0.000000
    0.003090  0.009511 0.000000
    0.000000  0.010000 0.000000
   -0.003090  0.009511 0.000000
   -0.005878  0.008090 0.000000
   -0.008090  0.005878 0.000000
   -0.009511  0.003090 0.000000
   -0.010000  0.000000 0.000000
   -0.009511 -0.003090 0.000000
   -0.008090 -0.005878 0.000000
   -0.005878 -0.008090 0.000000
   -0.003090 -0.009511 0.000000
   -0.000000 -0.010000 0.000000
    0.003090 -0.009511 0.000000
    0.005878 -0.008090 0.000000
    0.008090 -0.005878 0.000000
    0.009511 -0.003090 0.000000
}
```

`nband`: Denote the band number of which spin texture is calculated.

For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be three files in the `Out/Spin_Texture` folder, namely `kpt.dat` and `spin_texture.dat`, `plot_spin_texture.py`, corresponding to the k-point and the spin texture, the drawing script.
