# CHIRALITY

## introduction

Examines the chirality of Weyl points by calculating the Berry curvature on a sphere around the $\mathbf{k}$ point.

## example

An example is provided in the `examples/MnSb2Te4` directory to investigate the chirality of the Weyl point in MnSb$_2$Te$_4$.

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

CHIRALITY
{  
    k_vect                          0.0000 0.0000 -0.0538
    radius                          0.02  
    point_num                       100
}
```

`k_vect`: The k-point direct coordinates need to be calculated.

`radius`: The radius of the integrating sphere. The unit is $\AA^{-1}$ .

`point_num`: The number of k-points that are uniformly sampled on a spherical surface.

The `chirality.dat` file that record the results are in the `Out/Chirality` folder.
