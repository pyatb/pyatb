# PDOS

## introduction

The distribution of electronic states at various energies is characterized by the density of states (DOS), while the partial density of states (PDOS) is a useful tool for analyzing the contribution of individual atomic orbitals to the
DOS.

The implementation of PDOS is given by

$$
g_{\mu}(E)= \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \sum_{n} \sum_{\nu} C_{n\nu}^{*}(\mathbf{k}) S_{\nu\mu}(\mathbf{k}) C_{n\mu}(\mathbf{k}) \delta(E - E_{n\mathbf{k}}),
$$

where $C_{n\mu}(\mathbf{k})$ is the coefficient of the NAO.

## example

An example (refer to folder `examples/Si2`) of calculating the PDOS of the diamond Si is given here.

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

PDOS
{
    stru_file     STRU
    e_range       -5.0 17.0
    de            0.01
    sigma         0.07
    kpoint_mode   mp
    mp_grid       12 12 12
}
```

`stru_file`: The structure file name of the supercell. This file indicates the crystal structure of the supercell and the corresponding orbital file. Make sure that both the structure file and the orbital file exist.

`e_range`: Specify the energy range of dos, the unit is eV.

`de`: specifies the energy interval.

`sigma`: Specify the parameters of Gaussian smearing.

For the k point setting of this function, please refer to the `kpoint_mode` module.

After the task calculation is completed, there will be three files in the `Out/PDOS` folder, namely `TDOS.dat` and `PDOS.xml`, `plot_dos.py`. Specify the projected atomic orbital index in the plot script, and then draw the PDOS plot.
