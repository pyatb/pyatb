# BAND_STRUCTURE

An example of calculating the energy band structure of the topological insulator Bi$_2$Se$_3$ is provided in the folder `examples/Bi2Se3`.

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

BAND_STRUCTURE
{
    wf_collect                     0
    kpoint_mode                    line
    kpoint_num                     5
    high_symmetry_kpoint
    0.00000 0.00000 0.0000 100  # G
    0.00000 0.00000 0.5000 100  # Z
    0.50000 0.50000 0.0000 100  # F
    0.00000 0.00000 0.0000 100  # G
    0.50000 0.00000 0.0000 1    # L
}
```

`wf_collect`: Whether to output wave function matrix information. The wave function file stores the expansion coefficients of NAOs.

There are three ways to set k-points: k-point, k-line, and k-mesh, with the keyword `kpoint_mode` used to define the mode. The setting parameters for each mode differ, so please refer to the `INPUT` for detailed instructions. In this example, the `line` mode is used. In this mode, `kpoint_num` specifies the number of high-symmetry points, and `high_symmetry_kpoint` records the direct coordinates of these points and the number of k-points between each pair of high-symmetry points. Each row in the setting consists of four numbers, where the first three indicate the coordinates, and the last number specifies the number of k-points between the given k-point and the next high symmetry point.

After the calculation is done, a number of files will be generated in the `Out/Band_Structure` folder, including `kpt.dat` and `band.dat`. If `kpoint_mode` is set to `line`, the output will also include the energy band diagram and its corresponding drawing script.
