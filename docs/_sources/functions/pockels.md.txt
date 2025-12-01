# POCKELS

## introduction

Pockels effect which is a refractive index modulation created by the photoinduced space charge field combined with an applied static or nearly static external field.
The Pockels coefficient is defined:

$$
P^c(\omega) = \chi^{abc}(\omega_{\Sigma};\omega_1\rightarrow0,\omega_2) E^a(\omega_1) E^b(\omega_2)
$$



For whether the coupling between external field and light field is positive or negative is undetermined, we define in simplification:

$$
\chi(\omega_{\Sigma};\omega_1,\omega_2) = (\chi(\omega_1+\omega_2;\omega_1,\omega_2)+\chi(\omega_1-\omega_2;\omega_1,\omega_2))/2
$$



 In PYATB, the Pockels coefficient is calculated using the  form which divides Pockels into inter- and intra-band parts, shown in the following formula:

$$
\chi_{inter,c;a,b}^{(2)} = \frac{e^3}{2\hbar^2}  \sum_{\mathbf{k}, n m l} \frac{r_{n m}^c r_{m l}^a r_{l n}^b}{s_1 \omega_{l n}-s_2 \omega_{m l}}\left(\frac{f_{n m}}{\omega_s-\omega_{m n}}-\frac{s_2 f_{n l}}{\omega_2-\omega_{l n}}\right. 
         \left.+\frac{s_1 f_{m l}}{\omega_1-\omega_{m l}}\right)+\left(a, \omega_1 ; b, \omega_2\right)
$$

$$
\chi_{intra,c;a,b}^{(2)}=\frac{i e^3}{(\omega_1+\omega_2) \hbar^2} \sum_{\mathbf{k},n, m } \frac{v_{m n}^c}{\omega_{n m}-(\omega_1+\omega_2)}\left(\frac{r_{n m}^a f_{m n}}{\omega_{n m}-\omega_1}\right)_{; k^b}+\frac{v_{m n}^c}{\omega_{n m}-(\omega_1+\omega_2)}\left(\frac{r_{n m}^b f_{m n}}{\omega_{n m}-\omega_2}\right)_{; k^a}
$$



We employ this form guaranteeing the convergance of Pockels coefficient near zero frequency of external field. It should be noted that this formula only considers the electronic contribution of Pockels effect, yet the ignored ironic part and ele-phonon interaction part might be crucial in certain circumstances.

## example

Here, we provide an example of calculating the Pockels coefficient of the GaAs (refer to folder `examples/GaAs`).

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        7.6331377499
    fermi_energy_unit   eV
    HR_route            ../data-HR-sparse_SPIN0.csr
    SR_route            ../data-SR-sparse_SPIN0.csr
    rR_route            ../data-rR-sparse.csr
    HR_unit             Ry
    rR_unit             Bohr
   max_kpoint_num        100000
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
POCKELS
{
    omega1   0
    omega    0.01 6
    domega   0.01
    grid     30 30 30

}
```

`omega1`: To set the frequency of the external field in term of energy, you can adjust it. the unit is eV.

`omega`: To set the energy range for the light field, you can adjust it. the unit is eV.

`domega`: Specifies the energy interval of the omega.

`grid`: Specifies the uniform k-point grid used to calculate the SHG.

Once the task has been finished, two crucial files are produced in the `Out/Pockels` directory. These files consist of `pockels.dat` and `plot_pockels.py`. 
The first file stores the Pockels' magnitude. The second file contains the script used for generating the visualization of the Pockels coefficient.
