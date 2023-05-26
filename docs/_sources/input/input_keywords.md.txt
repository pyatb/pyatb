# Full List of INPUT Keywords

- [INPUT_PARAMETERS](#input_parameters)

  [nspin](#nspin-input_nspin) | [package](#package-input_package) | [fermi_energy](#fermi_energy-input_fermi_energy) | [fermi_energy_unit](#fermi_energy_unit-input_fermi_energy_unit) | [HR_route](#hr_route-input_hr_route) | [SR_route](#sr_route-input_sr_route) | [rR_route](#rr_route-input_rr_route) | [binary](#binary-input_binary) | [HR_unit](#hr_unit-input_hr_unit) | [rR_unit](#rr_unit-input_rr_unit) | [max_kpoint_num](#max_kpoint_num-input_max_kpoint_num) | [sparse_format](#sparse_format-input_sparse_format)

- [LATTICE](#lattice)

  [lattice_constant](#lattice_constant-lattice_lattice_constant) | [lattice_constant_unit](#lattice_constant_unit-lattice_lattice_constant_unit) | [lattice_vector](#lattice_vector-lattice_lattice_vector)

- [BAND_STRUCTURE](#band_structure)

  [wf_collect](#wf_collect-bandstructure_wf_collect) | [kpoint_mode](#kpoint_mode-bandstructure_kpoint_mode)

- [BANDUNFOLDING](#bandunfolding)

  [stru_file](#stru_file-bandunfolding_stru_file) | [ecut](#ecut-bandunfolding_ecut) | [band_range](#band_range-bandunfolding_band_range) | [m_matrix](#m_matrix-bandunfolding_m_matrix) | [kpoint_mode](#kpoint_mode-bandunfolding_kpoint_mode)

- [FERMI_ENERGY](#fermi_energy)

  [temperature](#temperature-fermienergy_temperature) | [electron_num](#electron_num-fermienergy_electron_num) | [grid](#grid-fermienergy_grid) | [epsilon](#epsilon-fermienergy_epsilon) | 

- [FERMI_SURFACE](#fermi_surface)

  [bar](#bar-fermisurface_bar)  | [nbands](#nbands-fermisurface_nbands) | [kpoint_mode](#kpoint_mode-fermisurface_kpoint_mode)

- [FIND_NODES](#find_nodes)

  [energy_range](#energy_range-findnodes_energy_range) | [k_start](#k_start-findnodes_k_start) | [k_vect1](#k_vect1-findnodes_k_vect1) | [k_vect2](#k_vect2-findnodes_k_vect2) | [k_vect3](#k_vect3-findnodes_k_vect3) | [initial_grid](#initial_grid-findnodes_initial_grid) | [initial_threshold](#initial_threshold-findnodes_initial_threshold) | [adaptive_grid](#adaptive_grid-findnodes_adaptive_grid)[adaptive_threshold](#adaptive_threshold-findnodes_adaptive_threshold) | [kpoint_mode](#kpoint_mode-findnodes_kpoint_mode)

- [PDOS](#pdos)

  [stru_file](#stru_file-pdos_stru_file) | [e_range](#e_range-pdos_e_range) | [de](#de-pdos_de) | [sigma](#sigma-pdos_sigma) | [kpoint_mode](#kpoint_mode-pdos_kpoint_mode)

- [FAT_BAND](#fat_band)

  [band_range](#band_range-fatband_band_range) | [stru_file](#stru_file-fatband_stru_file) | [kpoint_mode](#kpoint_mode-fatband_kpoint_mode)

- [SPIN_TEXTURE](#spin_texture)

  [nband](#nband-spintexture_nband) | [kpoint_mode](#kpoint_mode-spintexture_kpoint_mode)

- [WILSON_LOOP](#wilson_loop)

  [occ_band](#occ_band-wilsonloop_occ_band) | [k_start](#k_start-wilsonloop_k_start) | [k_vect1](#k_vect1-wilsonloop_k_vect1) | [k_vect2](#k_vect2-wilsonloop_k_vect2) | [nk1](#nk1-wilsonloop_nk1) | [nk2](#nk2-wilsonloop_nk2)

- [POLARIZATION](#polarization)

  [occ_band](#occ_band-polarization_occ_band) | [nk1](#nk1-polarization_nk1) | [nk2](#nk2-polarization_nk2) | [nk3](#nk3-polarization_nk3) | [atom_type](#atom_type-polarization_atom_type) | [stru_file](#stru_file-polarization_stru_file) | [valence_e](#valence_e-polarization_valence_e)

- [BERRY_CURVATURE](#berry_curvature)

  [method](#method-berrycurvature_method)  | [occ_band](#occ_band-berrycurvature_occ_band) | [kpoint_mode](#kpoint_mode-berrycurvature_kpoint_mode)

- [AHC](#ahc)

  [method](#method-ahc_method) | [integrate_mode](#integrate_mode-ahc_integrate_mode)

- [CHERN_NUMBER](#chern_number)

  [method](#method-chernnumber_method) | [occ_band](#occ_band-chernnumber_occ_band) | [k_start](#k_start-chernnumber_k_start) | [k_vect1](#k_vect1-chernnumber_k_vect1) | [k_vect2](#k_vect2-chernnumber_k_vect2) | [integrate_mode](#integrate_mode-chernnumber_integrate_mode)

- [CHIRALITY](#chirality)

  [method](#method-chirality_method) | [k_vect](#k_vect-chirality_k_vect) | [radius](#radius-chirality_radius) | [point_num](#point_num-chirality_point_num)

- [JDOS](#jdos)

  [occ_band](#occ_band-jdos_occ_band) | [omega](#omega-jdos_omega) | [domega](#domega-jdos_domega) | [eta](#eta-jdos_eta) | [grid](#grid-jdos_grid)

- [OPTICAL_CONDUCTIVITY](#optical_conductivity)

  [occ_band](#occ_band-opticalconductivity_occ_band) | [omega](#omega-opticalconductivity_omega) | [domega](#domega-opticalconductivity_domega) | [eta](#eta-opticalconductivity_eta) | [grid](#grid-opticalconductivity_grid)

- [SHIFT_CURRENT](#shift_current)

  [occ_band](#occ_band-shiftcurrent_occ_band) | [omega](#omega-shiftcurrent_omega) | [domega](#domega-shiftcurrent_domega) | [smearing_method](#smearing_method-shiftcurrent_smearing_method) | [eta](#eta-shiftcurrent_eta) | [grid](#grid-shiftcurrent_grid) | [method](#method-shiftcurrent_method)

- [BERRY_CURVATURE_DIPOLE](#berry_curvature_dipole)

  [omega](#omega-berrycurvaturedipole_omega) | [domega](#domega-berrycurvaturedipole_domega) | [integrate_mode](#integrate_mode-berrycurvaturedipole_integrate_mode)


[Setting of k points](#setting-of-k-points)

- [When kpoint_mode is 'mp'](#when-kpoint_mode-is-mp)

  [mp_grid](#mp_grid) | [k_start](#k_start) | [k_vect1](#k_vect1) | [k_vect2](#k_vect2) | [k_vect3](#k_vect3)

- [When kpoint_mode is 'line'](#when-kpoint_mode-is-line)

  [kpoint_num](#kpoint_num) | [high_symmetry_kpoint](#high_symmetry_kpoint)

- [When kpoint_mode is 'direct'](#when-kpoint_mode-is-direct)

  [kpoint_num](#kpoint_num-1) | [kpoint_direct_coor](#kpoint_direct_coor)

[Setting of integration](#setting-of-integration)

- [When integrate_mode is 'Grid'](#when-integrate_mode-is-grid)

  [integrate_grid](#integrate_grid) | [adaptive_grid](#adaptive_grid) | [adaptive_grid_threshold](#adaptive_grid_threshold)

- [When integrate_mode is 'Adaptive'](#when-integrate_mode-is-adaptive)

  [relative_error](#relative_error) | [absolute_error](#absolute_error) | [initial_grid](#initial_grid)

## INPUT_PARAMETERS

### nspin {#input_nspin}

- **Type**: Integer
- **Description**: Indicates the spin component of the wave function, related to the structure of the HR file.
  - 1: regardless of spin.
  - 2: the wave function is divided into two groups, one group is all up and one group is all down.
  - 4: the wave function has both up and down components.
- **Default**: No default value

### package {#input_package}

- **Type**: String
- **Description**: Indicates data sources for HR, SR, rR.
- **Default**: ABACUS

### fermi_energy {#input_fermi_energy}

- **Type**: Real
- **Description**: Indicates the Fermi energy of the system. When set to Auto, the FERMI_ENERGY function needs to be added.
- **Default**: Auto

### fermi_energy_unit {#input_fermi_energy_unit}

- **Type**: String
- **Description**: The unit of Fermi. Can be set to Ry, eV.
- **Default**: eV

### HR_route {#input_HR_route}

- **Type**: String
- **Description**: Path to HR matrix file. When nspin=2, two sets of paths need to be provided.
- **Default**: No default value

### SR_route {#input_SR_route}

- **Type**: String
- **Description**: Path to the SR matrix file.
- **Default**: No default value

### rR_route {#input_rR_route}

- **Type**: String
- **Description**: Path to the rR matrix file.
- **Default**: No default value
  
### binary {#input_binary}

- **Type**: Boolean
- **Description**: Whether HR, SR, and rR files are binary files.
- **Default**: 0

### HR_unit {#input_HR_unit}

- **Type**: String
- **Description**: The unit of HR. Can be set to Ry, eV.
- **Default**: Ry

### rR_unit {#input_rR_unit}

- **Type**: String
- **Description**: The unit of rR. Can be set to Bohr, Angstrom.
- **Default**: Bohr

### max_kpoint_num {#input_max_kpoint_num}

- **Type**: Integer
- **Description**: The upper limit of the number of k points stored in the memory during program calculation, which is used to control the memory consumption during calculation.
- **Default**: 8000

### sparse_format {#input_sparse_format}

- **Type**: Boolean
- **Description**: Whether HR, SR, rR matrices are stored in memory is sparse storage.
- **Default**: 0


## LATTICE

### lattice_constant {#lattice_lattice_constant}

- **Type**: Real
- **Description**: The lattice constant of the system.
- **Default**: No default value

### lattice_constant_unit {#lattice_lattice_constant_unit}

- **Type**: String
- **Description**: The unit of the lattice constant. Can be set to Bohr, Angstrom.
- **Default**: Bohr

### lattice_vector {#lattice_lattice_vector}

- **Type**: Real
- **Description**: The 3 lattice vectors of the system. Each lattice vector is a row, with a total of 3 rows and 9 parameters.
- **Default**: No default value

## BAND_STRUCTURE

### wf_collect {#bandstructure_wf_collect}

- **Type**: Boolean
- **Description**: Whether to output wave function matrix information.
- **Default**: No default value

### kpoint_mode {#bandstructure_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value

## BANDUNFOLDING

### stru_file {#bandunfolding_stru_file}

- **Type**: String
- **Description**: Specify the strucutre file path.
- **Default**: No default value

### ecut {#bandunfolding_ecut}

- **Type**: Real
- **Description**: Used to determine the number of plane wave basis sets. Unit is Ry.
- **Default**: 10

### band_range {#bandunfolding_band_range}

- **Type**: Integer
- **Description**: Specifies the range of supercell energy band index within which the energy bands will be calculated. There are two parameters, representing the starting band index and the end band index, the index counts from 1.
- **Default**: No default value

### m_matrix {#bandunfolding_m_matrix}

- **Type**: Real
- **Description**: The lattice vector transformation matrix between the supercell and the primitive cell, with 9 parameters, is written on the same line.
- **Default**: No default value

### kpoint_mode {#bandunfolding_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point of unitcell. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value


## FERMI_ENERGY

### temperature {#fermienergy_temperature}

- **Type**: Real
- **Description**: temperature. The unit is K
- **Default**: 0

### electron_num {#fermienergy_electron_num}

- **Type**: Integer
- **Description**: The total number of electrons in the system.
- **Default**: No default value

### grid {#fermienergy_grid}

- **Type**: Integer
- **Description**: The grid to use for Newton interpolation. There are three parameters.
- **Default**: 10 10 10

### epsilon {#fermienergy_epsilon}

- **Type**: Real
- **Description**: Newton interpolation parameters, absolute accuracy.
- **Default**: 0.001


## FERMI_SURFACE

### bar {#fermisurface_bar}

- **Type**: Real
- **Description**: The max tolerable error bar for the Fermi surface
- **Default**: No default value

### nbands {#fermisurface_nbands}

- **Type**: Integer
- **Description**: If you know the energy band range where the Fermi energy is located, you can set this parameter to speed up the calculation. There are two numbers in total, indicating the range of the energy band. The default value is `0 0`, that is, all energy bands are considered.
- **Default**: 0 0

### kpoint_mode {#fermisurface_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value


## FIND_NODES

### energy_range {#findnodes_energy_range}

- **Type**: Integer
- **Description**: The energy range in which the program searches for degenerate points, the energy unit is eV.
- **Default**: No default value

### k_start {#findnodes_k_start}

- **Type**: Real
- **Description**: The origin point coordinates used to describe a Brillouin zone plane.
- **Default**: 0.0 0.0 0.0
  
### k_vect1 {#findnodes_k_vect1}

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 1.0 0.0 0.0

### k_vect2 {#findnodes_k_vect2}

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 0.0 1.0 0.0

### k_vect3 {#findnodes_k_vect3}

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 0.0 0.0 1.0

### initial_grid {#findnodes_initial_grid}

- **Type**: Integer
- **Description**: Set the initial grid for searching degenerate k points. There are three parameters.
- **Default**: 10 10 10

### initial_threshold {#findnodes_initial_threshold}

- **Type**: Real
- **Description**: The energy unit is eV. In the initial grid, only the k-points whose band differences are less than the threshold can enter the search for the next round of degenerate points.
- **Default**: 0.1

### adaptive_grid {#findnodes_adaptive_grid}

- **Type**: Integer
- **Description**: The refined grid will refine the k-points that reach the initial_threshold in the initial grid. There are three parameters.
- **Default**: 20 20 20

### adaptive_threshold {#findnodes_adaptive_threshold}

- **Type**: Real
- **Description**: The minimum difference considered in independent bands, the energy unit is eV. This means if the band gap is below this bar, it will be recognized as degenerate bands.
- **Default**: 0.001

### kpoint_mode {#findnodes_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value


## PDOS

### stru_file {#pdos_stru_file}

- **Type**: String
- **Description**: The structure file name. This file records the structure of the lattice, the types of elements, and the atomic orbitals used. Make sure that both the structure file and the orbital file exist.
- **Default**: No default value

### e_range {#pdos_e_range}

- **Type**: Real
- **Description**: The range of energy E. There are two parameters, indicating the starting point and the ending point.
- **Default**: No default value

### de {#pdos_de}

- **Type**: Real
- **Description**: The interval dE for the energy E.
- **Default**: 0.01

### sigma {#pdos_sigma}

- **Type**: Real
- **Description**: Parameters for gauss smearing.
- **Default**: 0.001

### kpoint_mode {#pdos_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value


## FAT_BAND

### band_range {#fatband_band_range}

- **Type**: Integer
- **Description**: There are two numbers (separated by spaces) to indicate which bands are selected for projection, counting from 1.
- **Default**: No default value

### stru_file {#fatband_stru_file}

- **Type**: String
- **Description**: The structure file name. This file indicates the crystal structure and the corresponding orbital file. Make sure that both the structure file and the orbital file exist.
- **Default**: No default value

### kpoint_mode {#fatband_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point of unitcell. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value


## SPIN_TEXTURE

### nband {#spintexture_nband}

- **Type**: Integer
- **Description**: A band index. (Band index counts from 1)
- **Default**: No default value

### kpoint_mode {#spintexture_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value


## WILSON_LOOP

### occ_band {#wilsonloop_occ_band}

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator.
- **Default**: No default value

### k_start {#wilsonloop_k_start}

- **Type**: Real
- **Description**: The origin point coordinates used to describe a Brillouin zone plane.
- **Default**: 0.0 0.0 0.0
  
### k_vect1 {#wilsonloop_k_vect1}

- **Type**: Real
- **Description**: The expansion vector is a vector used to define a Brillouin zone plane, and it is also the direction of integration for calculations.
- **Default**: 1.0 0.0 0.0

### k_vect2 {#wilsonloop_k_vect2}

- **Type**: Real
- **Description**: The expansion vector is a vector used to define a Brillouin zone plane, and it is also the direction of Wilson loop evolution for calculations.
- **Default**: 0.0 1.0 0.0

### nk1 {#wilsonloop_nk1}

- **Type**: Integer
- **Description**: k_vect1 is divided into nk1 k-points.
- **Default**: 100

### nk2 {#wilsonloop_nk2}

- **Type**: Integer
- **Description**: k_vect2 is divided into nk2 k-points.
- **Default**: 100


## POLARIZATION

### occ_band {#polarization_occ_band}

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator.
- **Default**: No default value

### nk1 {#polarization_nk1}

- **Type**: Integer
- **Description**: The number of samples in the x direction of reciprocal lattice vector $\mathbf{G}$.
- **Default**: No default value

### nk2 {#polarization_nk2}

- **Type**: Integer
- **Description**: The number of samples in the y direction of reciprocal lattice vector $\mathbf{G}$.
- **Default**: No default value

### nk3 {#polarization_nk3}

- **Type**: Integer
- **Description**: The number of samples in the z direction of reciprocal lattice vector $\mathbf{G}$.
- **Default**: No default value

### atom_type {#polarization_atom_type}

- **Type**: Integer
- **Description**: The number of element types in the system.
- **Default**: No default value

### stru_file {#polarization_stru_file}

- **Type**: String
- **Description**: Specify the strucutre file. NAOs files are not required.
- **Default**: No default value

### valence_e {#polarization_valence_e}

- **Type**: Integer
- **Description**: The number of valence electrons per element.
- **Default**: No default value


## BERRY_CURVATURE

### method {#berrycurvature_method}

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

### occ_band {#berrycurvature_occ_band}

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator. When this value is not set, it will be determined according to the Fermi energy.
- **Default**: -1

### kpoint_mode {#berrycurvature_kpoint_mode}

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value


## AHC

### method {#ahc_method}

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

### integrate_mode {#ahc_integrate_mode}

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).
- **Default**: No default value


## CHERN_NUMBER

### method {#chernnumber_method}

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

### occ_band {#chernnumber_occ_band}

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator. When this value is not set, it will be determined according to the Fermi energy.
- **Default**: -1

### k_start {#chernnumber_k_start}

- **Type**: Real
- **Description**: The origin point coordinates used to describe a Brillouin zone plane.
- **Default**: 0.0 0.0 0.0
  
### k_vect1 {#chernnumber_k_vect1}

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 1.0 0.0 0.0

### k_vect2 {#chernnumber_k_vect2}

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 0.0 1.0 0.0

### integrate_mode {#chernnumber_integrate_mode}

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).
- **Default**: No default value

## CHIRALITY

### method {#chirality_method}

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

### k_vect {#chirality_k_vect}

- **Type**: Real
- **Description**: The k-point coordinates need to be calculated. There are three parameters to represent the coordinates.
- **Default**: No default value

### radius {#chirality_radius}

- **Type**: Real
- **Description**: The radius of the integrating sphere. The unit is $\AA^{-1}$ .
- **Default**: No default value

### point_num {#chirality_point_num}

- **Type**: Integer
- **Description**: The number of k-points that are uniformly sampled on a spherical surface.
- **Default**: No default value


## JDOS

### occ_band {#jdos_occ_band}

- **Type**: Integer
- **Description**: Specifies the occupied energy band of the system. Currently, only insulator or semiconductor materials can be calculated.
- **Default**: No default value

### omega {#jdos_omega}

- **Type**: Real
- **Description**: Specifies the photon energy, the unit is eV. There are two parameters, indicating the starting point and the ending point.
- **Default**: No default value

### domega {#jdos_domega}

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

### eta {#jdos_eta}

- **Type**: Real
- **Description**: Specify the parameters of Gaussian smearing.
- **Default**: 0.01

### grid {#jdos_grid}

- **Type**: Integer
- **Description**: The grid for integration. There are 3 parameters in total.
- **Default**: No default value


## OPTICAL_CONDUCTIVITY

### occ_band {#opticalconductivity_occ_band}

- **Type**: Integer
- **Description**: Used to specify the occupied energy band of an insulator or semiconductor. Currently this function can only calculate insulators or semiconductors.
- **Default**: No default value

### omega {#opticalconductivity_omega}

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

### domega {#opticalconductivity_domega}

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

### eta {#opticalconductivity_eta}

- **Type**: Real
- **Description**: Parameters for triangular smearing.
- **Default**: 0.01

### grid {#opticalconductivity_grid}

- **Type**: Integer
- **Description**: The grid for integration. There are 3 parameters in total.
- **Default**: No default value


## SHIFT_CURRENT

### occ_band {#shiftcurrent_occ_band}

- **Type**: Integer
- **Description**: Used to specify the occupied energy band of an insulator or semiconductor. Currently this function can only calculate insulators or semiconductors.
- **Default**: No default value

### omega {#shiftcurrent_omega}

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

### domega {#shiftcurrent_domega}

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

### smearing_method {#shiftcurrent_smearing_method}

- **Type**: Integer
- **Description**: The method of smearing. `0`: no smearing. `1`: Gaussian smearing. `2`: adaptive smearing.
- **Default**: 1

### eta {#shiftcurrent_eta}

- **Type**: Real
- **Description**: Specify the parameters of Gaussian smearing.
- **Default**: 0.01

### grid {#shiftcurrent_grid}

- **Type**: Integer
- **Description**: The grid for integration. There are 3 parameters in total.
- **Default**: No default value

### method {#shiftcurrent_method}

- **Type**: Integer
- **Description**: Specify the method to calculate the shift current. `0` represents calculation using the Sternheimer equation, `1` represents the first order partial derivative calculation.
- **Default**: 1


## BERRY_CURVATURE_DIPOLE

### omega {#berrycurvaturedipole_omega}

- **Type**: Real
- **Description**: To set the energy range for the Berry curvature dipole, you can adjust it based on the Fermi energy level. The unit is eV. There are two parameters.
- **Default**: No default value

### domega {#berrycurvaturedipole_domega}

- **Type**: Real
- **Description**: Specifies the energy interval of the omega.
- **Default**: No default value

### integrate_mode {#berrycurvaturedipole_integrate_mode}

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).Since the integration is of a tensor, only 'Grid'integrate_mode is available.
- **Default**: No default value


# Setting of k points

As long as the kpoint_mode parameter exists in FUNCTIONS, the following setting methods are to be followed.

## When kpoint_mode is 'mp'

### mp_grid

- **Type**: Integer
- **Description**: The grid dividing the Brillouin zone. There are three parameters to divide the three-dimensional Brillouin zone.
- **Default**: No default value

### k_start

- **Type**: Real
- **Description**: The origin point coordinates of the Brillouin zone.
- **Default**: 0.0 0.0 0.0

### k_vect1

- **Type**: Real
- **Description**: Expanded vector of the Brillouin zone.
- **Default**: 1.0 0.0 0.0

### k_vect2

- **Type**: Real
- **Description**: Expanded vector of the Brillouin zone.
- **Default**: 0.0 1.0 0.0
 
### k_vect3

- **Type**: Real
- **Description**: Expanded vector of the Brillouin zone.
- **Default**: 0.0 0.0 1.0

## When kpoint_mode is 'line'

### kpoint_num

- **Type**: Integer
- **Description**: The number of high symmetry points.
- **Default**: No default value


### high_symmetry_kpoint

- **Type**: Real
- **Description**: Fractional coordinates of high symmetry points and line densities of corresponding k-lines. The first three parameters are the fractional coordinates of the high symmetry points, and the fourth parameter is the line density.
- **Default**: No default value

## When kpoint_mode is 'direct'

### kpoint_num

- **Type**: Integer
- **Description**: the number of k points.
- **Default**: No default value

### kpoint_direct_coor

- **Type**: Real
- **Description**: Fractional coordinates of the k point.
- **Default**: No default value


# Setting of integration

As long as the integrate_mode parameter exists in FUNCTIONS, the following setting methods are followed.

## When integrate_mode is 'Grid'

### integrate_grid

- **Type**: Integer
- **Description**: Low precision grid for integration. There are three parameters.
- **Default**: 4 4 4

### adaptive_grid

- **Type**: Integer
- **Description**: High precision grid for integration. There are three parameters.
- **Default**: 4 4 4

### adaptive_grid_threshold

- **Type**: Real
- **Description**: If the value of a k point is greater than this value, then the k point will be adapted.
- **Default**: 50.0

## When integrate_mode is 'Adaptive'

### relative_error

- **Type**: Real
- **Description**: The relative error of the adaptive integral.
- **Default**: 1e-6

### absolute_error

- **Type**: Real
- **Description**: The absolute error of the adaptive integral.
- **Default**: 0.1

### initial_grid

- **Type**: Real
- **Description**: The initial grid for adaptive integration. There are three parameters.
- **Default**: 1 1 1
