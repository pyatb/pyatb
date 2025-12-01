# SPIN_TEXTURE

## introduction

Spin texture refers to the spatial distribution of electron spins in momentum space. In PYATB, the spin texture is calculated as follows,

$$
\langle \Psi_{n\mathbf{k}} | \hat{\sigma}_i  | \Psi_{n\mathbf{k}} \rangle
= \sum_{\mu,\nu,s,s^{\prime}} C^{*}_{n,\mu s}(\mathbf{k}) S_{\mu\nu, s s^{\prime}}(\mathbf{k}) \hat{\sigma}_{i,s s^{\prime}} C_{n,\nu s^{\prime}}(\mathbf{k}),
$$

where $\hat{\sigma}_i$ are the Pauli matrices, with $i$= $x$, $y$, $z$, and $s$=$\uparrow$, $\downarrow$ is the spin index.

## example

Here, we provide an example of calculating the spin texture of WSSe (refer to folder `tutorial/WSSe_spin_texture`). 

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                           4
    package                         ABACUS
    fermi_energy                    1.4232523137
    fermi_energy_unit               eV
    HR_route                        ../abacus/OUT.WSSe/data-HR-sparse_SPIN0.csr
    SR_route                        ../abacus/OUT.WSSe/data-SR-sparse_SPIN0.csr
    rR_route                        ../abacus/OUT.WSSe/data-rR-sparse.csr
    HR_unit                         Ry
    rR_unit                         Bohr
    max_kpoint_num                  4000
}

LATTICE
{
    lattice_constant                1.8897162
    lattice_constant_unit           Bohr
    lattice_vector
    3.2521424294         0.0000000000         0.0000000000
    1.6260712147         2.8164379605         0.0000000000
    0.0000000000         0.0000000000        18.2324867249
}

BAND_STRUCTURE
{
    wf_collect                0
    kpoint_mode               line
    kpoint_num                7
    high_symmetry_kpoint
    0.0000000000   0.0000000000   0.0000000000  100    # G
    0.0000000000   0.5000000000   0.0000000000  100    # M'
    0.3333333333   0.6666666666   0.0000000000  100    # K'
    0.0000000000   0.0000000000   0.0000000000  100    # G
    0.5000000000   0.5000000000   0.0000000000  100    # M
    0.6666666666   0.3333333333   0.0000000000  100    # K
    0.0000000000   0.0000000000   0.0000000000  1      # G
    kpoint_label                    G,M',K',G,M,K,G
}

SPIN_TEXTURE
{
    band_range                35 38
    kpoint_mode               line
    kpoint_num                7
    high_symmetry_kpoint
    0.0000000000   0.0000000000   0.0000000000  100    # G
    0.0000000000   0.5000000000   0.0000000000  100    # M'
    0.3333333333   0.6666666666   0.0000000000  100    # K'
    0.0000000000   0.0000000000   0.0000000000  100    # G
    0.5000000000   0.5000000000   0.0000000000  100    # M
    0.6666666666   0.3333333333   0.0000000000  100    # K
    0.0000000000   0.0000000000   0.0000000000  1      # G
    kpoint_label                    G,M',K',G,M,K,G
}
```

`band_range`: Integer type, with two values, used to define the range of energy bands for calculation. Counting starts from 1. For example, “1 20” specifies that the calculation includes energy bands from the 1st to the 20th, a total of 21 bands.

For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be three files in the `Out/Spin_Texture` folder, namely `kpt.dat` and `spin_texture_x.dat`, `spin_texture_y.dat`, `spin_texture_z.dat`, `plot_spintexture_line.py`, corresponding to the k-point and the spin texture, the drawing script.
