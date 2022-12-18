# Band structure

An example (refer to folder `example/Bi2Se3`) of calculating the energy band structure of the topological insulator Bi$_2$Se$_3$ is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin                          4
    package                        ABACUS
    fermi_energy                   9.481194886038594
    fermi_energy_unit              eV
    HR_route                       data-HR-sparse_SPIN0.csr
    SR_route                       data-SR-sparse_SPIN0.csr
    rR_route                       new-data-rR-tr_SPIN4
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

There are three ways to set kpoint, the keyword is `kpoint_mode` (refer to Inupt), after the calculation is done, several files will be generated, `kpt.dat`, `band.dat`. When `kpoint_mode` is set to `line`, the energy band diagram and the corresponding drawing script will be output.