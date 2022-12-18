# Find nodes

## introduction

Find nodes is a method of finding k points with degenerate energy bands in a given energy range.

## example

An example (refer to folder `example/MnBi2Te4`) of calculating the spectral function of the MnBi$_2$Te$_4$ slab is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin               4
    package             ABACUS
    fermi_energy        9.2284811591974805
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            new-data-rR-tr_SPIN4
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant   1
    lattice_constant_unit   Angstrom
    lattice_vector
    4.3773399999000002 0.0000000000000000 0.0000000000000000
    2.1886700000000001 3.7908876409999999 0.0000000000000000
    2.1886700000000001 1.2636292099999999 13.7730333333000008
}

FIND_NODES
{
       energy_range    9.1984811591974805   9.2584811591974805
       bar           1e-2
        kpoint_mode     mp
        k_start        -0.1 -0.1 -0.2
        k_vect1        0.2 0.0 0.0
        k_vect2        0.0 0.2 0.0
        k_vect3        0.0 0.0 0.4
        mp_grid        10 10 10
}
```
`energy_range`: The energy range in which the program searches for degenerate points, the energy unit is eV.

`bar`: The minimum difference considered in independent bands, the energy unit is eV. This means if the band gap is below this bar, it will be recognized as degenerate bands.

For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be three files in the `Out/Find_Nodes` folder, namely `nodes_kpt.dat` ,`plot_nodes.py`and`nodes.pdf`, corresponding to the degenerate k-point(s) in direct coordinate and a plotting script with its plot.
