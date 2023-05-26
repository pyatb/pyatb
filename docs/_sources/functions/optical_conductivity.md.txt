# OPTICAL_CONDUCTIVITY

## introduction

The frequency-dependent optical conductivity expressed by the Kubo-Greenwood formula can be formulated as

$$
\sigma_{\alpha\beta}(\hbar\omega) = -\frac{i e^2\hbar}{N V_{\mathrm{cell}}}\sum_{\mathbf{k}}
\sum_{n,m}\left(\frac{f_{n\mathbf{k}}-f_{m\mathbf{k}}}{E_{n\mathbf{k}}-E_{m\mathbf{k}}}\right)
\frac{\langle\psi_{n\mathbf{k}}|v_{\alpha}|\psi_{m\mathbf{k}}\rangle\langle\psi_{m\mathbf{k}}|v_{\beta}|\psi_{n\mathbf{k}}\rangle}{\hbar\omega + E_{n\mathbf{k}}-E_{m\mathbf{k}} + i\eta}\, .
$$


The imaginary part of the dielectric function is 

$$
\epsilon_i^{\alpha\beta}(\omega) = -\frac{e^2 \pi}{\epsilon_0 \hbar} \int \frac{d\mathbf{k}}{\left(2\pi\right)^3} \sum_{nm}f_{nm}r_{nm}^{\alpha}r_{mn}^{\beta} \delta\left(\omega_{mn} - \omega\right)\, ,
$$

The real part of the dielectric function is obtained by the Kramer-Kronig transformation,

$$
\epsilon_{r}^{\alpha\beta}(\omega) = \delta_{\alpha\beta} + \frac{2}{\pi} \mathbf{P} \int_{0}^{\infty} d\omega^{\prime} \frac{\omega^{\prime}\epsilon_{i}^{\alpha\beta}\left(\omega^{\prime}\right)}{\omega^{\prime 2} - \omega^2}\, .
$$


The linear optical spectrum can be calculated through the dielectric function, such as refractive index $n(\omega)$, extinction coefficient $\kappa(\omega)$, absorption coefficient $\alpha(\omega)$, energy-loss function $L(\omega)$, reflectivity $R(\omega)$:

$$
\begin{aligned}
n(\omega) &= \left[\frac{\sqrt{\varepsilon_1^2+\varepsilon_2^2}+\varepsilon_1}{2}\right]^{\frac{1}{2}} \\
\kappa(\omega) &= \left[\frac{\sqrt{\varepsilon_1^2+\varepsilon_2^2}-\varepsilon_1}{2}\right]^{\frac{1}{2}} \\
\alpha(\omega) &= \frac{\sqrt{2} \omega}{c}\left[\sqrt{\varepsilon_1^2+\varepsilon_2^2}-\varepsilon_1\right]^{\frac{1}{2}} \\
L(\omega) &= \operatorname{Im}\left(\frac{-1}{\varepsilon(\omega)}\right)=\frac{\varepsilon_2}{\varepsilon_1^2+\varepsilon_2^2} \\
R(\omega) &= \frac{(n-1)^2+k^2}{(n+1)^2+k^2}
\end{aligned}
$$

## example

Here, we provide an example (located in the `examples/Si` folder) demonstrating the calculation of the optical conductivity and dielectric function for diamond Si.

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        6.389728305291531
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            data-rR-sparse.csr
    HR_unit             Ry
    rR_unit             Bohr
}

LATTICE
{
    lattice_constant        1.8897162
    lattice_constant_unit   Bohr
    lattice_vector
    0.000000000000  2.715000000000  2.715000000000
    2.715000000000  0.000000000000  2.715000000000
    2.715000000000  2.715000000000  0.000000000000
}

OPTICAL_CONDUCTIVITY
{
    occ_band      4
    omega         0   10
    domega        0.01
    eta           0.1
    grid          50 50 50
}
```

`occ_band`: Used to specify the occupied energy band of an insulator or semiconductor. Currently this function can only calculate insulators or semiconductors.

`omega`: Specifies the photon energy, the unit is eV.

`domega`: Specifies the energy interval of the omega.

`eta`: Specify the parameters of Gaussian smearing.

`grid`: Specifies the uniform k-point grid used to calculate the optical conductivity.

After completing the task, five main files are generated in the `Out/Optical_Conductivity` folder, namely `optical_conductivity_real_part.dat`, `optical_conductivity_imag_part.dat`, `dielectric_function_real_part.dat`, 
`dielectric_function_imag_part.dat` and `plot_optical.py`.
