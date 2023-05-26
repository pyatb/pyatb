# FERMI_ENERGY

## introduction

This function is to calculate the Fermi energy of solid materials given temperature and electronic occupation number. If the system is an insulator, fermi energy is given by the valence band maximum (VBM) .

For each $\mathbf{k}$ point, the  probability of finding an electron in any energy state should obey the Fermi-Dirac distribution. The integration of occupied electrons over the entire Brillouin zone should be the occupation number. Though which, the exact Fermi energy could be attained following Newton interpolation.

$$
f(E,E_f,T)=\frac{1}{1+e^{\left(\frac{E-E_f}{k_B T}\right)}}
$$

$$
N[E_f]=\int_{BZ}[d\mathbf{k}]\sum_nf(E_n,E_f,T)
$$

## example

You can find an example of how to calculate the Fermi energy of Copper in the `examples/Cu` folder. 

The `Input` file is:

```
INPUT_PARAMETERS
{
    nspin                  1
    package                ABACUS
    fermi_energy           Auto
    fermi_energy_unit      eV
    HR_route               data-HR-sparse_SPIN0.csr
    SR_route               data-SR-sparse_SPIN0.csr
    HR_unit                Ry
    rR_unit                Bohr
    max_kpoint_num         8000
}

LATTICE
{
    lattice_constant       6.91640 
    lattice_constant_unit  Bohr
    lattice_vector
    0.50	0.50	0.00
    0.50	0.00	0.50
    0.00	0.50	0.50
}

FERMI_ENERGY
{
    temperature            0
    electron_num           11
    grid                   50 50 50
    epsilon                1e-4
}
```

`electron_num`: The number of the electrons in the system. The number of valence electrons, which is the sum of all atomic valence electrons of the system, can be obtained from the self-consistent output of the first-principles software.

`epsilon`: The max tolerable error of Newton interpolation. If two steps of Newton interpolation differs less than this epsilon, the calculation would stop and output the result.

`temperature`: The temperature of the system, unit in K.

After the task calculation is completed, there will be one file in the `Out/Fermi_Energy` folder, namely `fermi_energy.dat`, showing the calculated Fermi energy.
