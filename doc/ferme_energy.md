# Fermi energy

## introduction

Fermi energy is a function calculating the Fermi energy of a metal given temperature and electronic occupation number. Note that the Fermi energy of insulators might not be well-defined.

For each $\mathbf{k}$ point, the  probability of finding an electron in any energy state should obey the Fermi-Dirac distribution. The integration of occupied electrons over the entire Brillouin zone should be the occupation number. Though which, the exact Fermi energy could be attained following Newton interpolation.



$$
f(E,E_f,T)=\frac{1}{1+e^{\left(\frac{E-E_f}{k_B T}\right)}}
$$
$$
N[E_f]=\int_{BZ}[d\mathbf{k}]\sum_nf(E_n,E_f,T)
$$

## example

An example (refer to folder `example/Cu`) of calculating the spectral function of the Cu is given here.

The `Input` file are:

```txt {.line-numbers}
INPUT_PARAMETERS
{
    nspin               4
    package             ABACUS
    fermi_energy        9.540272009417667
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    HR_unit             Ry
    rR_unit             Bohr
    max_kpoint_num      8000
}

LATTICE
{
    lattice_constant        1.8897162
    lattice_constant_unit   Bohr
    lattice_vector
    4.3337001801 0 0 
    2.16685009 3.7530944483 0 
    4.3339143734 2.50218663 27.2734777344 
}

FERMI_ENERGY
{
	occupation_number
	temperature
    integration_mode   Grid
    integrate_grid     50 50 50 
    adaptive_grid      10 10 10 
    adaptive_threshold 1
    epsilon         1e-6
    
}
```
`occupation_number`: The number of the electrons in the system. It should be calculated by accumulating the orbital numbers  of each elements which is taken into consideration.

`epsilon`:The max tolerable error of Newton interpolation. If two steps of Newton interpolation differs less than this epsilon, the calculation would stop and output the answer.

`temperature`: The temperature of the system, unit in K.

After the task calculation is completed, there will be one file in the `Out/Fermi_Energy` folder, namely `fermi_energy.dat` ,showing the calculated Fermi energy. 
