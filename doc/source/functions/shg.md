# SECOND_HARMONIC_GENERATION

## introduction

Second harmonic generation (SHG), also called frequency doubling, is a nonlinear optical process, in which photons interacting with a nonlinear material are effectively ‘combined’ to form new photons having twice the frequency of initial photons.
The SHG coefficient is defined:

$$
P^c(2\omega) = \chi^{abc}(-2\omega;\omega,\omega) E^a(\omega) E^b(\omega)
$$



In PYATB, the SHG calculated using two methods. One is the commonly used form which divides SHG into inter- and intra-band parts, shown in the following formula:

$$
\chi_{interband}^{a b c}(-2 \omega, \omega, \omega)=  \frac{e^3}{\hbar^2 \Omega} \sum_{n m l, \mathbf{k}} \frac{r_{n m}^a\left\{r_{m l}^b r_{l n}^c\right\}}{\left(\omega_{l n}-\omega_{m l}\right)} 
\left[\frac{2 f_{n m}}{\omega_{m n}-2 \omega}+\frac{f_{l n}}{\omega_{l n}-\omega}+\frac{f_{m l}}{\omega_{m l}-\omega}\right]
$$

$$
\begin{aligned}
        \chi_{intraband}^{a b c}(-2 \omega, \omega, \omega)= & \frac{i}{2} \frac{e^3}{\hbar^2 \Omega} \sum_{n m, \mathbf{k}} f_{n m}\left[\frac{2}{\omega_{m n}\left(\omega_{m n}-2 \omega\right)} r_{n m}^a\left(r_{n m ; c}^b+r_{m n ; b}^c\right)+\frac{1}{\omega_{m n}\left(\omega_{m n}-\omega\right)}\left(r_{n m ; c}^a r_{m n}^b+r_{n m ; b}^a r_{m n}^c\right)\right. \\
        & +\frac{1}{\omega_{m n}^2}\left(\frac{1}{\omega_{m n}-\omega}-\frac{4}{\omega_{m n}-2 \omega}\right) r_{n m}^a\left(r_{m n}^b \Delta_{m n}^c+r_{m n}^c \Delta_{m n}^b\right)\\
        &-\left.\frac{1}{2 \omega_{m n}\left(\omega_{m n}-\omega\right)}\left(r_{n m ; a}^b r_{m n}^c+r_{n m ; a}^c r_{m n}^b\right)\right]
\end{aligned}
$$

This is the most widely used form and is the default method (method 0) of SHG calculation. Secondly SHG could be expressed in the form of multiplying velocity matrices:

$$
\chi^{abc}\left(-2\omega; \omega, \omega\right)=-\sum_{n, m, l,\mathbf{k}} \frac{i}{2\omega^3(2\omega-\omega_{m n})}\left(\frac{f_{n l}}{\omega-\omega_{l n}}+\frac{f_{m l}}{\omega-\omega_{m l}}\right)  v_{n m}^c v_{m l}^a v_{l n}^b
$$

This form is generally aligns with the upper formula, while the calculation can be costly for a wide energy range can be contributing to a single energy point. 

## example

Here, we provide an example of calculating the SHG of the GaAs (refer to folder `tutorial/GaAs_SHG/`).

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        10.171348972
    fermi_energy_unit   eV
    HR_route            ../abacus/OUT.GaAs/data-HR-sparse_SPIN0.csr
    SR_route            ../abacus/OUT.GaAs/data-SR-sparse_SPIN0.csr
    rR_route            ../abacus/OUT.GaAs/data-rR-sparse.csr
    HR_unit             Ry
    rR_unit             Bohr
   max_kpoint_num       100000
}

LATTICE
{
    lattice_constant        1.889727
    lattice_constant_unit   Bohr
    lattice_vector
    0.0000000000000000    2.7650000000000001    2.7650000000000001
    2.7650000000000001    0.0000000000000000    2.7650000000000001
    2.7650000000000001    2.7650000000000001    0.0000000000000000

}

SHG
{
    omega           0.01 4
    domega          0.01
    eta             0.05
    grid            50 50 50
}
```

`omega`: To set the energy range for the SHG, you can adjust it. the unit is eV.

`domega`: Specifies the energy interval of the omega.

`grid`: Specifies the uniform k-point grid used to calculate the SHG.

`eta`: $\hbar\omega \to \hbar\omega + i \eta$ is used to prevent numerical divergence caused by a zero denominator.

Once the task has been finished, three crucial files are produced in the `Out/Second_Harmonic_Generation` directory. These files consist of `shg_real.dat`, `shg_imag.dat` and `plot_shg.py`. 

The first two files contain the real and imaginary parts of the SHG's magnitude, and the last one is a plotting script.

