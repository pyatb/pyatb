# Input file <!-- omit in toc -->

- [Structure of the Input file](#structure-of-the-input-file)
- [List of keywords](#list-of-keywords)
  - [INPUT_PARAMETERS](#input_parameters)
    - [nspin](#nspin)
    - [package](#package)
    - [fermi_energy](#fermi_energy)
    - [fermi_energy_unit](#fermi_energy_unit)
    - [HR_route](#hr_route)
    - [SR_route](#sr_route)
    - [rR_route](#rr_route)
    - [binary](#binary)
    - [HR_unit](#hr_unit)
    - [rR_unit](#rr_unit)
    - [max_kpoint_num](#max_kpoint_num)
    - [sparse_format](#sparse_format)
  - [LATTICE](#lattice)
    - [lattice_constant](#lattice_constant)
    - [lattice_constant_unit](#lattice_constant_unit)
    - [lattice_vector](#lattice_vector)
  - [BAND_STRUCTURE](#band_structure)
    - [wf_collect](#wf_collect)
    - [kpoint_mode](#kpoint_mode)
  - [BANDUNFOLDING](#bandunfolding)
    - [stru_file](#stru_file)
    - [ecut](#ecut)
    - [band_range](#band_range)
    - [m_matrix](#m_matrix)
    - [kpoint_mode](#kpoint_mode-1)
  - [FERMI_ENERGY](#fermi_energy-1)
    - [temperature](#temperature)
    - [electron_num](#electron_num)
    - [grid](#grid)
    - [epsilon](#epsilon)
  - [FERMI_SURFACE](#fermi_surface)
    - [bar](#bar)
    - [nbands](#nbands)
    - [kpoint_mode](#kpoint_mode-2)
  - [FIND_NODES](#find_nodes)
    - [energy_range](#energy_range)
    - [initial_grid](#initial_grid)
    - [initial_threshold](#initial_threshold)
    - [adaptive_grid](#adaptive_grid)
    - [adaptive_threshold](#adaptive_threshold)
    - [k_start](#k_start)
    - [k_vect1](#k_vect1)
    - [k_vect2](#k_vect2)
    - [k_vect3](#k_vect3)
  - [JDOS](#jdos)
    - [occ_band](#occ_band)
    - [omega](#omega)
    - [domega](#domega)
    - [eta](#eta)
    - [grid](#grid-1)
  - [PDOS](#pdos)
    - [e_range](#e_range)
    - [de](#de)
    - [sigma](#sigma)
    - [kpoint_mode](#kpoint_mode-3)
    - [SPIN_TEXTURE](#spin_texture)
    - [nband](#nband)
    - [kpoint_mode](#kpoint_mode-4)
  - [AHC](#ahc)
    - [method](#method)
    - [integrate_mode](#integrate_mode)
  - [BERRY_CURVATURE](#berry_curvature)
    - [method](#method-1)
    - [occ_band](#occ_band-1)
    - [kpoint_mode](#kpoint_mode-5)
  - [BERRY_CURVATURE_DIPOLE](#berry_curvature_dipole)
    - [omega](#omega-1)
    - [domega](#domega-1)
    - [integrate_mode](#integrate_mode-1)
  - [CHERN_NUMBER](#chern_number)
    - [method](#method-3)
    - [occ_band](#occ_band-2)
    - [k_start](#k_start-2)
    - [k_vect1](#k_vect1-2)
    - [k_vect2](#k_vect2-2)
    - [integrate_mode](#integrate_mode-2)
  - [CHIRALITY](#chirality)
    - [method](#method-4)
    - [k_vect](#k_vect)
    - [radius](#radius)
    - [point_num](#point_num)
  - [CPGE](#cpge)
    - [omega](#omega-2)
    - [domega](#domega-2)
    - [integrate_mode](#integrate_mode-3)
  - [OPTICAL_CONDUCTIVITY](#optical_conductivity)
    - [occ_band](#occ_band-3)
    - [omega](#omega-3)
    - [domega](#domega-3)
    - [eta](#eta-1)
    - [grid](#grid-2)
  - [POLARIZATION](#polarization)
    - [occ_band](#occ_band-4)
    - [nk1](#nk1)
    - [nk2](#nk2)
    - [nk3](#nk3)
    - [atom_type](#atom_type)
    - [stru_file](#stru_file-1)
    - [valence_e](#valence_e)
  - [SHIFT_CURRENT](#shift_current)
    - [occ_band](#occ_band-5)
    - [omega](#omega-4)
    - [domega](#domega-4)
    - [eta](#eta-2)
    - [grid](#grid-3)
    - [method](#method-5)
  - [WILSON_LOOP](#wilson_loop)
    - [occ_band](#occ_band-6)
    - [k_start](#k_start-3)
    - [k_vect1](#k_vect1-3)
    - [k_vect2](#k_vect2-3)
    - [nk1](#nk1-1)
    - [nk2](#nk2-1)
  - [DRUDE_WEIGHT](#drude_weight)
    - [omega](#omega-5)
    - [domega](#domega-5)
    - [integrate_mode](#integrate_mode-4)
- [Setting of k points](#setting-of-k-points)
  - [1. When kpoint_mode is mp](#1-when-kpoint_mode-is-mp)
    - [mp_grid](#mp_grid)
    - [k_start](#k_start-4)
    - [k_vect1](#k_vect1-4)
    - [k_vect2](#k_vect2-4)
    - [k_vect3](#k_vect3-2)
  - [2. When kpoint_mode is line](#2-when-kpoint_mode-is-line)
    - [kpoint_num](#kpoint_num)
    - [high_symmetry_kpoint](#high_symmetry_kpoint)
  - [3. When kpoint_mode is direct](#3-when-kpoint_mode-is-direct)
    - [kpoint_num](#kpoint_num-1)
    - [kpoint_direct_coor](#kpoint_direct_coor)
- [Setting of integration](#setting-of-integration)
  - [1. When integrate_mode is Grid](#1-when-integrate_mode-is-grid)
    - [integrate_grid](#integrate_grid)
    - [adaptive_grid](#adaptive_grid-1)
    - [adaptive_grid_threshold](#adaptive_grid_threshold)
  - [2. When integrate_mode is Adaptive](#2-when-integrate_mode-is-adaptive)
    - [relative_error](#relative_error)
    - [absolute_error](#absolute_error)
    - [initial_grid](#initial_grid-1)

## Structure of the Input file

In this file, all text field starting with `#` or `\\` are considered comments. All parameters are lowercase. The Input file is mainly decomposed into three parts: `INPUT_PARAMETERS`, `LATTICE`, `FUNCTIONS`. The `INPUT_PARAMETERS` part mainly describes the basic information of the system and the environment settings of the execution program, `LATTICE` describes the structure information of the system, and `FUNCTIONS` is the parameter setting of each function itself. The three parts are composed as follows:

```text {.line-numbers}
INPUT_PARAMETERS
{
    ...  // Parameters
}

LATTICE
{
    ... // Parameters
}

FUNCTIONS
{
    ... // Parameters
}
```

A simple Input file is as follows:

```text {.line-numbers}
INPUT_PARAMETERS
{
    nspin               4
    package             ABACUS
    fermi_energy        9.540272009417667
    fermi_energy_unit   eV
    HR_route            data-HR-sparse_SPIN0.csr
    SR_route            data-SR-sparse_SPIN0.csr
    rR_route            new-data-rR-tr_SPIN4
    HR_unit             Ry
    rR_unit             Bohr
    max_kpoint_num      8000
}

LATTICE
{
    lattice_constant        1.8897162
    lattice_constant_unit   Bohr
    lattice_vector
    -2.069  -3.583614  0.000000
     2.069  -3.583614  0.000000
     0.000   2.389075  9.546667
}

BAND_STRUCTURE  // FUNCTIONS 1
{
    nbands                  1  98
    wf_collect              1
    kpoint_mode             line
    kpoint_num              5
    high_symmetry_kpoint
    0.00000 0.00000 0.0000 100
    0.00000 0.00000 0.5000 100
    0.50000 0.50000 0.0000 100
    0.00000 0.00000 0.0000 100
    0.50000 0.00000 0.0000 1
}

PDOS            // FUNCTIONS 2
{
    e_range                 3 15
    de                      0.01
    sigma                   0.1
    kpoint_mode             mp
    mp_grid                 10 10 10
}
```


## List of keywords

### INPUT_PARAMETERS

#### nspin

- **Type**: Integer
- **Description**: Indicates the spin component of the wave function, related to the structure of the HR file.
  - 1: regardless of spin.
  - 2: the wave function is divided into two groups, one group is all up and one group is all down.
  - 4: the wave function has both up and down components.
- **Default**: No default value

#### package

- **Type**: String
- **Description**: Indicates data sources for HR, SR, rR.
- **Default**: ABACUS

#### fermi_energy

- **Type**: Real
- **Description**: Indicates the Fermi energy of the system. When set to Auto, the FERMI_ENERGY function needs to be added.
- **Default**: Auto

#### fermi_energy_unit

- **Type**: String
- **Description**: The unit of Fermi. Can be set to Ry, eV.
- **Default**: eV

#### HR_route

- **Type**: String
- **Description**: Path to HR matrix file. When nspin=2, two sets of paths need to be provided.
- **Default**: No default value

#### SR_route

- **Type**: String
- **Description**: Path to the SR matrix file.
- **Default**: No default value

#### rR_route

- **Type**: String
- **Description**: Path to the rR matrix file.
- **Default**: No default value
  
#### binary

- **Type**: Boolean
- **Description**: Whether HR, SR, and rR files are binary files.
- **Default**: 0

#### HR_unit

- **Type**: String
- **Description**: The unit of HR. Can be set to Ry, eV.
- **Default**: Ry

#### rR_unit

- **Type**: String
- **Description**: The unit of rR. Can be set to Bohr, Angstrom.
- **Default**: Bohr

#### max_kpoint_num

- **Type**: Integer
- **Description**: The upper limit of the number of k points stored in the memory during program calculation, which is used to control the memory consumption during calculation.
- **Default**: 8000

#### sparse_format

- **Type**: Boolean
- **Description**: Whether HR, SR, rR matrices are stored in memory is sparse storage.
- **Default**: 0


### LATTICE

#### lattice_constant

- **Type**: Real
- **Description**: The lattice constant of the system.
- **Default**: No default value

#### lattice_constant_unit

- **Type**: String
- **Description**: The unit of the lattice constant. Can be set to Bohr, Angstrom.
- **Default**: Bohr

#### lattice_vector

- **Type**: Real
- **Description**: The 3 lattice vectors of the system. Each lattice vector is a row, with a total of 3 rows and 9 parameters.
- **Default**: No default value

### BAND_STRUCTURE

#### wf_collect

- **Type**: Boolean
- **Description**: Whether to output wave function matrix information.
- **Default**: No default value

#### kpoint_mode

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value

### BANDUNFOLDING

#### stru_file

- **Type**: String
- **Description**: Specify the strucutre file path.
- **Default**: No default value

#### ecut

- **Type**: Real
- **Description**: Used to determine the number of plane wave basis sets. Unit is eV.
- **Default**: 10

#### band_range

- **Type**: Integer
- **Description**: Specifies the range of supercell energy band index within which the energy bands will be calculated. There are two parameters, representing the starting band index and the end band index, the index counts from 1.
- **Default**: No default value

#### m_matrix

- **Type**: Real
- **Description**: The lattice vector transformation matrix between the supercell and the primitive cell, with 9 parameters, is written on the same line.
- **Default**: No default value

#### kpoint_mode

- **Type**: String
- **Description**: Used to set the k point of unitcell. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value

### FERMI_ENERGY

#### temperature

- **Type**: Real
- **Description**: temperature. The unit is K
- **Default**: 0

#### electron_num

- **Type**: Integer
- **Description**: The total number of electrons in the system.
- **Default**: No default value

#### grid

- **Type**: Integer
- **Description**: The grid to use for Newton interpolation. There are three parameters.
- **Default**: 10 10 10

#### epsilon

- **Type**: Real
- **Description**: Newton interpolation parameters, absolute accuracy.
- **Default**: 0.001

### FERMI_SURFACE

#### bar

- **Type**: 
- **Description**: 
- **Default**: 

#### nbands

- **Type**: 
- **Description**: 
- **Default**: 

#### kpoint_mode

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value

### FIND_NODES

#### energy_range

- **Type**: 
- **Description**: 
- **Default**: 

#### initial_grid

- **Type**: 
- **Description**: 
- **Default**: 

#### initial_threshold

- **Type**: 
- **Description**: 
- **Default**: 

#### adaptive_grid

- **Type**: 
- **Description**: 
- **Default**: 

#### adaptive_threshold

- **Type**: 
- **Description**: 
- **Default**: 

#### k_start

- **Type**: 
- **Description**: 
- **Default**: 

#### k_vect1

- **Type**: 
- **Description**: 
- **Default**: 

#### k_vect2

- **Type**: 
- **Description**: 
- **Default**: 

#### k_vect3

- **Type**: 
- **Description**: 
- **Default**: 

### JDOS

#### occ_band

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator.
- **Default**: No default value

#### omega

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

#### domega

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

#### eta

- **Type**: Real
- **Description**: Parameters for gauss smearing.
- **Default**: 0.01

#### grid

- **Type**: Integer
- **Description**: The grid for integration. There are 3 parameters in total.
- **Default**: No default value

### PDOS

#### e_range

- **Type**: Real
- **Description**: The range of energy E. There are two parameters, indicating the starting point and the ending point.
- **Default**: No default value

#### de

- **Type**: Real
- **Description**: The interval dE for the energy E.
- **Default**: 0.01

#### sigma

- **Type**: Real
- **Description**: Parameters for gauss smearing.
- **Default**: 0.001

#### kpoint_mode

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value

#### SPIN_TEXTURE

#### nband

- **Type**: Integer
- **Description**: A band index. (Band index counts from 1)
- **Default**: No default value

#### kpoint_mode

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value

### AHC

#### method

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

#### integrate_mode

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).
- **Default**: No default value

### BERRY_CURVATURE

#### method

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

#### occ_band

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator. When this value is not set, it will be determined according to the Fermi energy.
- **Default**: -1

#### kpoint_mode

- **Type**: String
- **Description**: Used to set the k point. See [Setting of k points](#setting-of-k-points)
- **Default**: No default value

### BERRY_CURVATURE_DIPOLE
#### omega

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

#### domega

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

#### integrate_mode

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).Since the integration is of a tensor, only 'Grid'integrate_mode is available.
- **Default**: No default value


### CHERN_NUMBER

#### method

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

#### occ_band

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator. When this value is not set, it will be determined according to the Fermi energy.
- **Default**: -1

#### k_start

- **Type**: Real
- **Description**: The origin point coordinates used to describe a Brillouin zone plane.
- **Default**: 0.0 0.0 0.0
  
#### k_vect1

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 1.0 0.0 0.0

#### k_vect2

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 0.0 1.0 0.0

#### integrate_mode

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).
- **Default**: No default value

### CHIRALITY

#### method

- **Type**: Integer
- **Description**: Method for calculating berry curvature. `0` means direct calculation, `1` means calculation by Kubo formula.
- **Default**: 0

#### k_vect

- **Type**: Real
- **Description**: The k-point coordinates need to be calculated. There are three parameters to represent the coordinates.
- **Default**: No default value

#### radius

- **Type**: Real
- **Description**: Indicates the radius.
- **Default**: No default value

#### point_num

- **Type**: Integer
- **Description**: Calculate the number of k-points required for chirality.
- **Default**: No default value

### CPGE

#### omega

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

#### domega

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

#### integrate_mode

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).Since the integration is of a tensor, only 'Grid'integrate_mode is available.
- **Default**: No default value


### OPTICAL_CONDUCTIVITY

#### occ_band

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator.
- **Default**: No default value

#### omega

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

#### domega

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

#### eta

- **Type**: Real
- **Description**: Parameters for triangular smearing.
- **Default**: 0.01

#### grid

- **Type**: Integer
- **Description**: The grid for integration. There are 3 parameters in total.
- **Default**: No default value

### POLARIZATION

#### occ_band

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator.
- **Default**: No default value

#### nk1

- **Type**: Integer
- **Description**: The number of samples in the x direction of reciprocal lattice vector $\mathbf{G}$.
- **Default**: No default value

#### nk2

- **Type**: Integer
- **Description**: The number of samples in the y direction of reciprocal lattice vector $\mathbf{G}$.
- **Default**: No default value

#### nk3

- **Type**: Integer
- **Description**: The number of samples in the z direction of reciprocal lattice vector $\mathbf{G}$.
- **Default**: No default value

#### atom_type

- **Type**: Integer
- **Description**: types of elements in the system.
- **Default**: No default value

#### stru_file

- **Type**: String
- **Description**: Specify the strucutre file path.
- **Default**: No default value

#### valence_e

- **Type**: Integer
- **Description**: The number of valence electrons per element.
- **Default**: No default value

### SHIFT_CURRENT

#### occ_band

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator.
- **Default**: No default value

#### omega

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

#### domega

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

#### eta

- **Type**: Real
- **Description**: Parameters for gauss smearing.
- **Default**: 0.01

#### grid

- **Type**: Integer
- **Description**: The grid for integration. There are 3 parameters in total.
- **Default**: No default value

#### method

- **Type**: Integer
- **Description**: Specify the method to calculate the shift current. `0` represents calculation using the Sternheimer equation, `1` represents the first order partial derivative calculation.
- **Default**: 1

### WILSON_LOOP

#### occ_band

- **Type**: Integer
- **Description**: The number of occupied energy bands of an insulator.
- **Default**: No default value

#### k_start

- **Type**: Real
- **Description**: The origin point coordinates used to describe a Brillouin zone plane.
- **Default**: 0.0 0.0 0.0
  
#### k_vect1

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 1.0 0.0 0.0

#### k_vect2

- **Type**: Real
- **Description**: The expansion vector used to describe a Brillouin zone plane.
- **Default**: 0.0 1.0 0.0

#### nk1

- **Type**: Integer
- **Description**: k_vect1 is divided into nk1 k-points.
- **Default**: 100

#### nk2

- **Type**: Integer
- **Description**: k_vect2 is divided into nk2 k-points.
- **Default**: 100

### Drude Weight

#### omega

- **Type**: Real
- **Description**: The range of $\omega$. There are two parameters, indicating the starting point and the ending point. Unit is eV.
- **Default**: No default value

#### domega

- **Type**: Real
- **Description**: The energy interval of $\omega$.
- **Default**: No default value

#### integrate_mode

- **Type**: String
- **Description**: Used for integration settings. See [Setting of integration](#setting-of-integration).Since the integration is of a tensor, only 'Grid'integrate_mode is available.
- **Default**: No default value

## Setting of k points

As long as the kpoint_mode parameter exists in FUNCTIONS, the following setting methods are to be followed.

### 1. When kpoint_mode is 'mp'

#### mp_grid

- **Type**: Integer
- **Description**: The grid dividing the Brillouin zone. There are three parameters to divide the three-dimensional Brillouin zone.
- **Default**: No default value

#### k_start

- **Type**: Real
- **Description**: The origin point coordinates of the Brillouin zone.
- **Default**: 0.0 0.0 0.0

#### k_vect1

- **Type**: Real
- **Description**: Expanded vector of the Brillouin zone.
- **Default**: 1.0 0.0 0.0

#### k_vect2

- **Type**: Real
- **Description**: Expanded vector of the Brillouin zone.
- **Default**: 0.0 1.0 0.0
 
#### k_vect3

- **Type**: Real
- **Description**: Expanded vector of the Brillouin zone.
- **Default**: 0.0 0.0 1.0

### 2. When kpoint_mode is 'line'

#### kpoint_num

- **Type**: Integer
- **Description**: The number of high symmetry points.
- **Default**: No default value


#### high_symmetry_kpoint

- **Type**: Real
- **Description**: Fractional coordinates of high symmetry points and line densities of corresponding k-lines. The first three parameters are the fractional coordinates of the high symmetry points, and the fourth parameter is the line density.
- **Default**: No default value

### 3. When kpoint_mode is 'direct'

#### kpoint_num

- **Type**: Integer
- **Description**: the number of k points.
- **Default**: No default value

#### kpoint_direct_coor

- **Type**: Real
- **Description**: Fractional coordinates of the k point.
- **Default**: No default value


## Setting of integration

As long as the integrate_mode parameter exists in FUNCTIONS, the following setting methods are followed.

### 1. When integrate_mode is Grid

#### integrate_grid

- **Type**: Integer
- **Description**: Low precision grid for integration. There are three parameters.
- **Default**: 4 4 4

#### adaptive_grid

- **Type**: Integer
- **Description**: High precision grid for integration. There are three parameters.
- **Default**: 4 4 4

#### adaptive_grid_threshold

- **Type**: Real
- **Description**: If the value of a k point is greater than this value, then the k point will be adapted.
- **Default**: 50.0

### 2. When integrate_mode is Adaptive

#### relative_error

- **Type**: Real
- **Description**: The relative error of the adaptive integral.
- **Default**: 1e-6

#### absolute_error

- **Type**: Real
- **Description**: The absolute error of the adaptive integral.
- **Default**: 0.1

#### initial_grid

- **Type**: Real
- **Description**: The initial grid for adaptive integration. There are three parameters.
- **Default**: 1 1 1
