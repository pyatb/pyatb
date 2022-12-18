# Optical conductivity

# introduction

The frequency-dependent optical conductivity expressed by the Kubo-Greenwood formula can be formulated as

$$
\sigma_{\alpha\beta}(\hbar\omega) = -\frac{i e^2\hbar}{N V_{\mathrm{cell}}}\sum_{\mathbf{k}}
\sum_{n,m}\left(\frac{f_{n\mathbf{k}}-f_{m\mathbf{k}}}{E_{n\mathbf{k}}-E_{m\mathbf{k}}}\right)
\frac{\langle\psi_{n\mathbf{k}}|v_{\alpha}|\psi_{m\mathbf{k}}\rangle\langle\psi_{m\mathbf{k}}|v_{\beta}|\psi_{n\mathbf{k}}\rangle}{\hbar\omega + E_{n\mathbf{k}}-E_{m\mathbf{k}} + i\eta}
$$

# example

An example (refer to folder `example/Si`) of calculating the optical conductivity of the diamond Si is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin               1
    package             ABACUS
    fermi_energy        7.0
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            new-data-rR-tr_SPIN1
    HR_unit             Ry
    rR_unit             Bohr
    }

LATTICE
{
    lattice_constant        1.8897
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
    eta           0.2
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