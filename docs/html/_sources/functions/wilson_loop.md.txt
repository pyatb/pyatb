# WILSON_LOOP

## introduction

We can determine the Z2 topology number of the topological insulator by Wilson loop. Computing the six time reversal invariant planes, we can obtain the Z2 topology metrics ($\nu_0$, $\nu_1$ $\nu_2$ $\nu_3$). These six planes are k1=0.0, k1=0.5, k2=0.0, k2=0.5, k3=0.0, k3=0.5, respectively. When selecting a plane, only half of the plane needs to be selected.

$$
\begin{aligned}
\nu_0 &= Z2(ki=0) + Z2(ki=0.5) \quad mod \quad 2 \\
\nu_i &= Z2(ki=0.5)
\end{aligned}
$$

where $i = 1, 2, 3$ refers to the x, y and z directions.

The Wilson loop is implemented as follows:
$$
W_n(\mathbf{k_2}) = \frac{i}{2\pi} \int_{0}^{2\pi} d\mathbf{k_1} \langle u_{n,\mathbf{k_1}, \mathbf{k_2}} | \partial_{\mathbf{k_1}} | u_{n,\mathbf{k_1}, \mathbf{k_2}}\rangle.
$$


## example

Here is an example of how to calculate the Wilson loop of the topological insulator Bi$_2$Se$_3$ located in the `examples/Bi2Se3` folder.

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

WILSON_LOOP
{
    occ_band           78
    k_start            0.0  0.0  0.5
    k_vect1            1.0  0.0  0.0
    k_vect2            0.0  0.5  0.0
    nk1                101
    nk2                101
}
```

`occ_band`: The number of occupied energy bands of an insulator.

`k_start`: The origin point coordinates used to describe a Brillouin zone plane.

`k_vect1`: The expansion vector is a vector used to define a Brillouin zone plane, and it is also the direction of integration for calculations.

`k_vect2`: The expansion vector is a vector used to define a Brillouin zone plane, and it is also the direction of Wilson loop evolution for calculations. 

`nk1`: k_vect1 is divided into nk1 k-points.

`nk2`: k_vect2 is divided into nk2 k-points.

After the task calculation is completed, there will be two files in the `Out/Wilson_Loop` folder, namely `wilson_loop.dat` and `plot_wl.py`, corresponding to the Wilson loop of each k-point in the `k_vect2` direction and drawing script.
