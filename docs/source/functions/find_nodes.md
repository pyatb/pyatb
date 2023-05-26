# FIND_NODES

## introduction

Find nodes is a method of finding k points with degenerate energy bands in a given energy range. This function can be used to find the Weyl/Dirac points in Weyl/Dirac semimetals.

## example

An example of how to locate the Weyl point of MnSb$_2$Te$_4$ is provided in the `examples/MnSb2Te4` folder.

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                           4   
    package                         ABACUS
    fermi_energy                    9.9700823666762375
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
     4.2885731815169859   -0.0001203117360283   -0.0000592245216563
     2.1441823977387200    3.7140734770500492   -0.0000592245216895
     2.1439936829776851    1.2378353300160352   13.4147375436875418
}

FIND_NODES
{
    energy_range                    9.870 10.070
    k_start                         0.0 0.0 -0.2
    k_vect1                         0.0 0.0  0.0
    k_vect2                         0.0 0.0  0.0
    k_vect3                         0.0 0.0  0.4
    initial_grid                    1  1  100
    initial_threshold               0.01
    adaptive_grid                   1  1  20
    adaptive_threshold              0.001
}
```

`energy_range`: The energy range in which the program searches for degenerate points, the energy unit is eV.

Set the search space of k points: Selecting a parallel hexahedron in k-space requires an origin (`k_start`) and three vectors (`k_vect1`, `k_vect2`, `k_vect3`) that are not parallel to each other in general. In this example, `k_vect2` and `k_vect3` are zero vectors, so the chosen search space is k-line from (0.0, 0.0, -0.2) to (0.0, 0.0, 0.4).

After the task calculation is completed, there will be three files in the `Out/Find_Nodes` folder, namely `nodes.dat`, `plot_nodes.py` and `nodes.pdf`, corresponding to the degenerate k-point(s) in direct coordinate and a plotting script with its plot.
