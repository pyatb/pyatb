# SHIFT_CURRENT

## introduction

The shift current is an intrinsic contribution to the bulk photovoltaic effect (BPVE). It describes the photocurrent generated by light illumination on homogeneous non-centrosymmetric crystals. The shift current is a second-order optical response. It can be expressed as a DC current, generated by a monochromatic photoelectric field $\mathbf{E}(t) = \mathbf{E}(\omega)\mathrm{e}^{i\omega t} + \mathbf{E}(-\omega)\mathrm{e}^{-i\omega t}$, where

$$
J^{a} = 2 \sigma^{abc}(0; \omega, -\omega) E_{b}(\omega) E_{c}(-\omega).
$$

Here, $a, b, c = x, y, z$, and $\sigma^{abc}(0; \omega, -\omega)$ is the shift current tensor,

$$
\sigma^{abc}(0 ; \omega,-\omega) =	\frac{\pi e^3}{\hbar^2} \int \frac{d \mathbf{k}}{8 \pi^3} \sum_{n, m} f_{nm} \mathrm{Im}\left[ I_{mn}^{abc} + I_{mn}^{acb} \right] \delta\left(\omega_{m n}-\omega\right)
$$

where $I_{mn}^{abc} = r_{mn}^{b} r_{nm;a}^{c}$, $r_{nm}^a$ is the inter-band dipole matrix, and $r_{nm; a}^b$ is the generalized derivative of the dipole matrix.

Detailed descriptions can be found in [Second-order optical response in semiconductors](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.61.5337).

## example

Here, we provide an example of calculating the shift current conductivity of the monolayer WS$_2$ (refer to folder `examples/WS2`). 

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                   1
    package                 ABACUS
    fermi_energy            0.3484302262859574
    fermi_energy_unit       eV
    HR_route                data-HR-sparse_SPIN0.csr
    SR_route                data-SR-sparse_SPIN0.csr
    rR_route                data-rR-sparse.csr
    HR_unit                 Ry
    rR_unit                 Bohr
}

LATTICE
{
    lattice_constant        1.8897162
    lattice_constant_unit   Bohr
    lattice_vector
    3.183820900165   0.0              0.0            
   -1.591910450082   2.757269780643   0.0            
    0.0              0.0              20.086904001384
}

SHIFT_CURRENT
{
    occ_band                13
    omega                   0   4
    domega                  0.01
    smearing_method         1
    eta                     0.1
    grid                    1000 1000 1
}
```

`occ_band`: Used to specify the occupied energy band of an insulator or semiconductor. Currently this function can only calculate insulators or semiconductors.

`omega`: Specifies the photon energy, the unit is eV.

`domega`: Specifies the energy interval of the omega.

`eta`: Specify the parameters of Gaussian smearing.

`grid`: Specifies the uniform k-point grid used to calculate the shift current.

Upon completion of the task, two essential files are generated in the `Out/Shift_Current` folder. These files include `shift_current.dat` and `plot_shift_current.py`.
