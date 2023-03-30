# BERRY_CURVATURE_DIPOLE

## introduction

In a system with time-reversal symmetry, the Berry curvature is an odd function of $\mathbf{k}$, i.e., $\Omega_{a}(\mathbf{k})=-\Omega_{a}(-\mathbf{k})$. As a result, the integration of the Berry curvature over the BZ is zero. However, if the system lacks a inversion symmetry, a higher-order nonlinear AHC can arise. More specifically, $j_a^0=\chi_{abc} E_b (\omega) E_c(-\omega)$ and $j_a^{2 \omega}=\chi_{abc} E_b (\omega)E_c(\omega)$, describe a rectified current and the second harmonic, respectively, whereas $\omega$ is the driving frequency. The coefficient $\chi_{abc}$ is given by

$$
\chi_{abc}=-\varepsilon_{a d c} \frac{e^3 \tau}{2(1+i \omega \tau)} D_{bd} .
$$

where 

$$
D_{ab}=\int_{k} f_{0}\left(\frac{\partial \Omega_{b}}{\partial k_{a}} \right)
$$

is called the Berry curvature dipole. 

In PYATB, the Berry curvature dipole at a given temperature $T$ and chemical potential $\mu$ is calculated using the following formula:

$$
D_{a b}(\mu, T) = -\int \frac{\partial f_{0}(E,\mu, T)}{\partial E} D_{a b}(E) dE .
$$

## example

Here, we provide an example of calculating the Berry curvature dipole of the trigonal Te (refer to folder `examples/Te`).

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                     4
    package                   ABACUS
    fermi_energy              9.574476774876963
    fermi_energy_unit         eV
    HR_route                  data-HR-sparse_SPIN0.csr
    SR_route                  data-SR-sparse_SPIN0.csr
    rR_route                  data-rR-sparse.csr
    HR_unit                   Ry
    rR_unit                   Bohr
    max_kpoint_num            28000
}

LATTICE
{
    lattice_constant          1.8897162
    lattice_constant_unit     Bohr
    lattice_vector
    2.22    -3.84515    0.000
    2.22     3.84515    0.000
    0.00     0.00000    5.910

}

BERRY_CURVATURE_DIPOLE
{
    omega                     9.474  10.074
    domega                    0.001
    integrate_mode            Grid
    integrate_grid            400 400 400
    adaptive_grid             20 20 20
    adaptive_grid_threshold   20000
}
```

`omega`: To set the energy range for the Berry curvature dipole, you can adjust it based on the Fermi energy level. the unit is eV.

`domega`: Specifies the energy interval of the omega.

For k-point integration, please refer to the [Setting of integration] section.

Once the task has been finished, three crucial files are produced in the `Out/Berry_Curvature_Dipole` directory. These files consist of `bcd_step2.dat`, `kpoint_list`, and `plot_bcd.py`. The first file stores the Berry curvature dipole's magnitude, while the second file records the refined k-point coordinates. The third file contains the script used for generating the visualization of the dipole.
