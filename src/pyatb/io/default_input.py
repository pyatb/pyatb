"""
If the parameter is None, it indicates that this parameter has no default value and must be set.
If the parameter is [], it indicates that this parameter has no default value and don't need to be set in some conditions.

Each parameter definition follows the pattern:
    [python_type, expected_length_or_shape, default_value]

Examples:
- [int, 1, 0] means a scalar integer with default value 0
- [float, 3, [0.0, 0.0, 0.0]] means a 3-vector of floats with a default value
- [float, 3, 3, None] means a 3x3 float matrix without a default value
"""

# Function switches indicate which top-level blocks are available in an Input file.
function_switch = {
    # Fundamental and Shared Components for All Modules
    'INPUT_PARAMETERS'                     : False,
    'LATTICE'                              : False,

    # Band modules
    'BAND_STRUCTURE'                       : False,
    'BANDUNFOLDING'                        : False,
    'BANDUNFOLDING_SPIN_TEXTURE'           : False,
    'COHP'                                 : False,
    'FAT_BAND'                             : False,
    'FERMI_ENERGY'                         : False,
    'FERMI_SURFACE'                        : False,
    'FIND_NODES'                           : False,
    'PDOS'                                 : False,
    'SPIN_TEXTURE'                         : False,
    'SURFACE_STATE'                        : False,

    # Geometry modules
    'AHC'                                  : False,
    'ANC'                                  : False,
    'BERRY_CURVATURE'                      : False,
    'CHERN_NUMBER'                         : False,
    'CHIRALITY'                            : False,
    'ORBITAL_MAGNETIZATION'                : False,
    'POLARIZATION'                         : False,
    'SHC'                                  : False,
    'WILSON_LOOP'                          : False,
    
    # Optical modules
    'BERRY_CURVATURE_DIPOLE'               : False,
    'CPGE'                                 : False,
    'DRUDE_WEIGHT'                         : False,
    'JDOS'                                 : False,
    'OPTICAL_CONDUCTIVITY'                 : False,
    'POCKELS'                              : False,
    'SHG'                                  : False,
    'SHIFT_CURRENT'                        : False,

    # transport modules
    'BOLTZ_TRANSPORT'                      : False,
    
    # other utilities
    'REDUCE_BASIS'                         : False,
    
}

# Blocks listed here can be omitted even when empty because every field has a safe default.
block_can_be_empty = []

# These function blocks require the real-space position matrix `rR_route`.
need_rR_matrix = [
    'AHC',
    'ANC',
    'BERRY_CURVATURE',
    'BERRY_CURVATURE_DIPOLE',
    'CHERN_NUMBER',
    'CHIRALITY',
    'OPTICAL_CONDUCTIVITY',
    'ORBITAL_MAGNETIZATION',
    'POLARIZATION',
    'SHIFT_CURRENT',
    'WILSON_LOOP',
    'CPGE',
    'DRUDE_WEIGHT',
    'BOLTZ_TRANSPORT',
    'SHC',
    'SHG',
    'POCKELS'
]

# Sub-options attached to `kpoint_mode`.
kpoint_mode = {
    'mp' : 
    {
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],  # Origin of the sampled Brillouin-zone region
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],  # First spanning vector of the sampled k-space region
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],  # Second spanning vector of the sampled k-space region
        'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],  # Third spanning vector of the sampled k-space region
        'mp_grid'                     : [int, 3, None]                # Monkhorst-Pack-like divisions along the three spanning vectors
    },

    'line': 
    {
        'kpoint_num'                  : [int, 1, None],               # Number of high-symmetry points defining the path
        'high_symmetry_kpoint'        : [float, None, 4, None],       # Fractional coordinates plus line density for each path node
      # 'kpoint_num_in_line'          : [int, None, None]
        'kpoint_label'                : [str, None, None]             # Labels shown for the high-symmetry points on the plot
    },

    'direct' : 
    {
        'kpoint_num'                  : [int, 1, None],               # Number of explicitly listed k points
        'kpoint_direct_coor'          : [float, None, 3, None]        # Fractional coordinates of the explicit k points
    }
}

# Sub-options attached to `cal_surface_method` in `SURFACE_STATE`.
cal_surface_method = {
    'direct_diag' : 
    {
        'slab_layers'                 : [int, 1, None]                # Number of slab layers used in direct diagonalization
    },

    'direct_green' : 
    {
        'slab_layers'                 : [int, 1, None],               # Number of slab layers used to build the surface Green function
        'green_eta'                   : [float, 1, 0.001]             # Broadening in the Green-function denominator
    },

    'green_fun' : 
    {
        'green_eta'                   : [float, 1, 0.001],            # Broadening in the iterative surface Green-function method
    }
}

# Sub-options attached to `integrate_mode` in Berry and transport modules.
integrate_mode = {
    'Grid' : 
    {
        'integrate_grid'              : [int, 3, [4, 4, 4]],  # Coarse integration grid
        'adaptive_grid'               : [int, 3, [4, 4, 4]],  # Refined integration grid used near large-valued points
        'adaptive_grid_threshold'     : [float, 1, 50.0]      # Points above this value trigger adaptive refinement
    },

    'Adaptive' :
    {
        'relative_error'              : [float, 1, 1e-6],     # Target relative tolerance for adaptive integration
        'absolute_error'              : [float, 1, 0.1],      # Target absolute tolerance for adaptive integration
        'initial_grid'                : [int, 3, [1, 1, 1]]   # Starting grid before adaptive refinement
    }
}

# Source-package-specific route and unit settings.
package = {
    'ABACUS' :
    {
        'HR_route'                    : [str, None, None],    # Hamiltonian sparse matrix path(s); multiplicity depends on `nspin`
        'SR_route'                    : [str, 1, None],       # Overlap sparse matrix path for ABACUS-based workflows
        'rR_route'                    : [str, 1, []],         # Position-matrix path; only needed for modules in `need_rR_matrix`
        'binary'                      : [int, 1, 0],          # Whether matrix files are stored in binary format
        'HR_unit'                     : [str, 1, 'Ry'],       # Unit of the Hamiltonian matrix elements
        'rR_unit'                     : [str, 1, 'Bohr'],     # Unit of the position matrix elements
    },

    'WANNIER90' :
    {
        'w90_TB_route'                : [str, 1, None],       # Path to the Wannier90 tight-binding data file
        'w90_TB_has_r'                : [int, 1, 0]           # Whether the Wannier90 TB data already contains r-matrix information
    }
}

INPUT = {
    # Fundamental and Shared Components for All Modules
    'INPUT_PARAMETERS' : 
    {
        'nspin'                       : [int, 1, None],       # Spin mode of the imported model; also controls `HR_route` multiplicity
        'package'                     : [str, 1, 'ABACUS'],   # Source package for the TB data, e.g. ABACUS or WANNIER90
        'fermi_energy'                : [str, 1, 'Auto'],     # Fermi level or auto-detection flag used by downstream modules
        'fermi_energy_unit'           : [str, 1, 'eV'],       # Unit associated with `fermi_energy`
        'max_kpoint_num'              : [int, 1, 8000],       # Internal upper bound for handled k points
        'sparse_format'               : [int, 1, 0]           # Sparse-matrix storage/reading mode
    },

    'LATTICE' : 
    {
        'lattice_constant'            : [float, 1, None],     # Global lattice scaling factor
        'lattice_constant_unit'       : [str, 1, 'Bohr'],     # Unit of `lattice_constant`
        'lattice_vector'              : [float, 3, 3, None]   # 3x3 lattice-vector matrix
    },

    # Band modules
    'BAND_STRUCTURE' : 
    {
        'wf_collect'                  : [int, 1, False],      # Whether to collect wavefunction-related data for postprocessing
        'band_range'                  : [int, 2, [-1, -1]],   # Selected band window; [-1, -1] means use the default full range
        'kpoint_mode'                 : [str, 1, None]        # Path or list specification; usually `line` for band plots
    },

    'BANDUNFOLDING' : 
    {
        'stru_file'                   : [str, 1, None],       # Structure file used to reconstruct orbital/projector information
        'ecut'                        : [float, 1, 10],       # Plane-wave cutoff-like control used in the unfolding formalism
        'band_range'                  : [int, 2, None],       # Band window to unfold
        'm_matrix'                    : [float, 9, None],     # 3x3 unfolding transformation matrix written as 9 numbers
        'kpoint_mode'                 : [str, 1, None]        # Usually `line` for unfolded band paths
    },

    'BANDUNFOLDING_SPIN_TEXTURE' : 
    {
        'stru_file'                   : [str, 1, None],       # Structure file used for unfolding and spin projection metadata
        'ecut'                        : [float, 1, 10],       # Plane-wave cutoff-like control in the unfolding formalism
        'band_range'                  : [int, 2, None],       # Supercell band window to unfold and analyze
        'm_matrix'                    : [float, 9, None],     # 3x3 transformation matrix between supercell and primitive cell
        'kpoint_mode'                 : [str, 1, None]        # Unit-cell k-point path or list
    },

    'FAT_BAND' : 
    {
        'band_range'                  : [int, 2, None],       # Band window for orbital projection
        'stru_file'                   : [str, 1, None],       # Structure file required for orbital labels and projections
        'kpoint_mode'                 : [str, 1, None]        # Usually `line`
    },

    'COHP' :
    {
        'stru_file'                   : [str, 1, None],       # ABACUS STRU file used to resolve atoms and orbitals
        'atom_i_index'                : [int, 1, -1],         # First atom index in the ABACUS atom order, counted from 1
        'atom_j_index'                : [int, 1, -1],         # Second atom index in the ABACUS atom order, counted from 1
        'atom_i_orbs'                 : [str, -1, 'all'],     # Orbital selector for atom i
        'atom_j_orbs'                 : [str, -1, 'all'],     # Orbital selector for atom j
        'method'                      : [str, 1, 'COHP'],     # COHP or COOP population method
        'spin'                        : [str, 1, 'sum'],      # sum, up, or down
        'e_range'                     : [float, 2, None],     # Energy window for the projected population spectrum
        'de'                          : [float, 1, 0.05],     # Energy-grid spacing in eV
        'sigma'                       : [float, 1, 0.15],     # Gaussian smearing width in eV
        'invert'                      : [int, 1, 1],          # Whether to invert the COHP sign convention in output
        'shift_to_efermi'             : [int, 1, 1],          # Whether to shift energies by the Fermi level
        'output_prefix'               : [str, 1, 'COHP'],     # Prefix for COHP output files
        'input_file'                  : [str, 1, ''],         # Original Input path used for reproducibility metadata
        'orbital_dir'                 : [str, 1, ''],         # Directory containing numerical orbital files
        'kpoint_mode'                 : [str, 1, None],       # K-point specification for COHP sampling
        'kpoint_weights'              : [float, -1, []]     # Optional explicit k-point weights, normalized internally
    },

    'FERMI_ENERGY' : 
    {
        'temperature'                 : [float, 1, 0.0],           # Electronic temperature in K for Fermi-level determination
        'electron_num'                : [int, 1, None],            # Total number of electrons in the system
        'grid'                        : [int, 3, [10, 10, 10]],    # k-grid used by the Newton interpolation search
        'epsilon'                     : [float, 1, 1e-3]           # Absolute convergence threshold for the interpolation
    },

    'FERMI_SURFACE' : 
    {
        'bar'                         : [float, 1, 1e-3],    # Maximum tolerable error bar when locating the Fermi surface
        'nbands'                      : [int, 2, [0,0]],     # Optional band-index range near the Fermi level to accelerate the search
        'kpoint_mode'                 : [str, 1, None]       # Usually `mp` for 3D Fermi-surface sampling
    },

    'FIND_NODES':
    {
        'energy_range'                : [float, 2, [0.0, 0.0]],       # Target energy window for locating candidate crossings
        'initial_grid'                : [int, 3, [10, 10, 10]],       # Coarse search mesh
        'initial_threshold'           : [float, 1, 0.1],              # Energy-gap threshold on the coarse search
        'adaptive_grid'               : [int, 3, [20, 20, 20]],       # Refined mesh near candidate nodes
        'adaptive_threshold'          : [float, 1 , 1e-3],            # Tighter threshold for adaptive refinement
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],  # Origin of the searched k-space box
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],  # First spanning vector of the search box
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],  # Second spanning vector of the search box
        'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]]   # Third spanning vector of the search box
    },

    'PDOS' : 
    {
        'stru_file'                   : [str, 1, None],       # Structure file used to map orbitals and atoms
        'e_range'                     : [float, 2, None],     # Energy window for the DOS spectrum
        'de'                          : [float, 1, 0.01],     # Energy-grid spacing
        'sigma'                       : [float, 1, 1e-3],     # Broadening width for DOS smoothing
        'kpoint_mode'                 : [str, 1, None]        # K-point specification for the DOS sampling
    },

    'SPIN_TEXTURE' : 
    {
        'band_range'                  : [int, 2, None],       # Band window for reported spin expectation values
        'kpoint_mode'                 : [str, 1, None],       # Usually `line` or `direct`
        'orbital_resolved'            : [int, 1, 0],          # Whether to decompose the spin texture into orbital channels
        'stru_file'                   : [str, 1, 'STRU'],     # Structure metadata used for orbital-resolved output
    },

    'SURFACE_STATE':
    {
        'cal_surface_method'          : [str, 1, 'green_fun'],   # Surface-state solver choice
        'surface_direction'           : [str, 1, 'c'],           # Surface normal chosen from lattice-vector directions a, b, or c
        'energy_windows'              : [float, 2, [-1.0, 1.0]], # Energy window centered on the reference Fermi energy
        'de'                          : [float, 1, 0.01],        # Energy step inside the surface-state window
        'eta'                         : [float, 1, 0.01],        # Green-function broadening in eV
        'coupling_layers'             : [int, 1, None],          # Number of principal layers ensuring only nearest-layer coupling
        'calculate_layer'             : [int, 1, 1],             # Number of layers for which the spectral function is output
        'kpoint_mode'                 : [str, 1, None]           # Surface Brillouin-zone path or grid
    },

    # Geometry modules
    'AHC' : 
    {
        'method'                      : [int, 1, 0],          # Berry-curvature evaluation method: direct or Kubo-style
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_mode'              : [str, 1, None]        # Integration strategy in the Brillouin zone
    },

    'ANC' : 
    {
        'fermi_range'                 : [float, 2, [-1.0, 1.0]], # Fermi-level scan window relative to the reference Fermi energy
        'de'                          : [float, 1, 0.01],        # Energy step used in the Fermi-level scan
        'eta'                         : [float, 1, 0.01],        # Broadening in eV
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_grid'              : [int, 3, [4, 4, 4]]      # Brillouin-zone integration grid
    },

    'BERRY_CURVATURE' : 
    {
        'method'                      : [int, 1, 0],          # Berry-curvature evaluation algorithm selector
        'occ_band'                    : [int, 1, -1],         # Occupied-band cutoff; -1 usually means infer from the model or Fermi level
        'kpoint_mode'                 : [str, 1, None]        # Sampling mode for the requested curvature output
    },

    'CHERN_NUMBER' : 
    {
        'method'                      : [int, 1, 0],                  # Berry-curvature evaluation method: direct or Kubo-style
        'occ_band'                    : [int, 1, -1],                 # Occupied-band cutoff; infer from Fermi level when unset
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],  # Origin of the integration plane
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],  # First in-plane vector
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],  # Second in-plane vector
        'integrate_mode'              : [str, 1, None]                # Integration strategy for the Chern-number evaluation
    },

    'CHIRALITY':
    {
        'method'                      : [int, 1, 0],                 # Chirality-evaluation method selector
        'occ_band'                    : [int, 1, -1],                # Occupied-band cutoff for the enclosing surface calculation
        'k_vect'                      : [float, 3, [0.0, 0.0, 0.0]], # Node-centered direction or coordinate selector
        'radius'                      : [float, 1, 0.01],            # Radius of the enclosing sphere in reciprocal space
        'point_num'                   : [int, 1, 1000]               # Number of sampled points on the enclosing sphere
    },

    'ORBITAL_MAGNETIZATION' : 
    {
        'fermi_energy'                : [float, 1, None],        # Reference Fermi energy used for the orbital-magnetization scan
        'fermi_range'                 : [float, 2, [-2.0, 2.0]], # Energy window around the reference Fermi energy
        'de'                          : [float, 1, 0.05],        # Energy step in the scan
        'eta'                         : [float, 1, 0.01],        # Broadening in eV
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None],          # Brillouin-zone integration grid
    },

    'POLARIZATION' : 
    {
        'occ_band'                    : [int, 1, None],      # Number of occupied bands in the insulating ground state
        'nk1'                         : [int, 1, 8],         # Sampling count along reciprocal-lattice direction 1
        'nk2'                         : [int, 1, 8],         # Sampling count along reciprocal-lattice direction 2
        'nk3'                         : [int, 1, 8],         # Sampling count along reciprocal-lattice direction 3
        'atom_type'                   : [int, 1, None],      # Number of distinct atomic species in the structure
        'stru_file'                   : [str, 1, None],      # Structure file containing species and valence information
        'valence_e'                   : [int, None, None]    # Number of valence electrons for each atomic species
    },

    'SHC' : 
    {
        'alpha'                       : [str, 1, "x"],           # Charge-current direction in sigma_{alpha beta}^{gamma}
        'beta'                        : [str, 1, "y"],           # Electric-field direction in sigma_{alpha beta}^{gamma}
        'gamma'                       : [str, 1, "z"],           # Spin-polarization direction in sigma_{alpha beta}^{gamma}
        'fermi_range'                 : [float, 2, [-1.0, 1.0]], # Fermi-level scan window relative to the reference Fermi energy
        'de'                          : [float, 1, 0.01],        # Energy step in the Fermi-level scan
        'eta'                         : [float, 1, 0.01],        # Kubo broadening in eV
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_grid'              : [int, 3, [4, 4, 4]]      # Brillouin-zone integration grid
    },

    'WILSON_LOOP' : 
    {
        'occ_band'                    : [int, 1, None],               # Number of occupied bands included in the Wilson loop
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],  # Origin of the Brillouin-zone plane
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],  # Integration direction of the Wilson loop
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],  # Evolution direction of the Wilson loop
        'nk1'                         : [int, 1, 100],                # Number of points along the Wilson-loop direction
        'nk2'                         : [int, 1, 100],                # Number of base points sweeping the manifold
    },

    # Optical modules
    'BERRY_CURVATURE_DIPOLE':
    {
        'omega'                       : [float, 2, None],    # Energy window around the Fermi level used in the dipole scan
        'domega'                      : [float, 1, None],    # Energy step in the dipole scan
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None]       # Brillouin-zone integration grid
    },

    'CPGE':
    {
        'omega'                       : [float, 2, None],    # Photon-energy window for the circular photogalvanic response
        'domega'                      : [float, 1, None],    # Photon-energy step
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_mode'              : [str, 1, None]       # Integration strategy over the Brillouin zone
    },

    'DRUDE_WEIGHT':
    {
        'omega'                       : [float, 2, None],    # Frequency window for the Drude-weight spectrum
        'domega'                      : [float, 1, None],    # Frequency step
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_mode'              : [str, 1, None]       # Integration strategy over the Brillouin zone
    },

    'JDOS' : 
    {
        'occ_band'                    : [int, 1, None],      # Occupied-band cutoff for interband-transition counting
        'omega'                       : [float, 2, None],    # Photon-energy window in eV
        'domega'                      : [float, 1, None],    # Photon-energy step
        'eta'                         : [float, 1, 0.01],    # Gaussian-smearing width in eV
        'grid'                        : [int, 3, None]       # Brillouin-zone integration grid
    },

    'OPTICAL_CONDUCTIVITY' : 
    {
        'occ_band'                    : [int, 1, -1],           # Occupied-band cutoff used in the optical transition sum
        'omega'                       : [float, 2, [0.0, 4.0]], # Frequency window
        'domega'                      : [float, 1, 0.01],       # Frequency step
        'eta'                         : [float, 1, 0.01],       # Broadening in eV
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None],         # Dense BZ grid for response integration
        'method'                      : [int, 1, 1],            # Optical-conductivity evaluation method selector
        'static_dielectric_only'      : [int, 1, 0]             # If enabled, only the static dielectric response is computed
    },

    'POCKELS':
    {
        'omega1'                      : [float, 1, 0],       # Frequency of the external electric field in eV
        'omega'                       : [float, 2, None],    # Frequency window of the optical response in eV
        'domega'                      : [float, 1, None],    # Frequency step
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None]       # Brillouin-zone integration grid
    },

    'SHG':
    {
        'method'                      : [int, 1, 0],          # SHG evaluation method selector
        'eta'                         : [float, 1, 0.05],     # Broadening in eV
        'omega'                       : [float, 2, None],     # Frequency window
        'domega'                      : [float, 1, None],     # Frequency step
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None]        # Dense k-grid for the SHG response
    },

    'SHIFT_CURRENT' : 
    {
        'occ_band'                    : [int, 1, None],       # Occupied-band cutoff for nonlinear optical transitions
        'omega'                       : [float, 2, None],     # Frequency window
        'domega'                      : [float, 1, None],     # Frequency step
        'smearing_method'             : [int, 1, 1],          # 0: none, 1: Gaussian, 2: adaptive
        'eta'                         : [float, 1, 0.01],     # Broadening in eV
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None],       # Dense BZ grid
        'method'                      : [int, 1, 1],          # Shift-current evaluation method
        'n_occ'                       : [int, 1, -1],         # Occupied-band index selector, counted from 1
        'm_unocc'                     : [int, 1, -1]          # Unoccupied-band index selector, counted from 1
    },

    # transport modules
    'BOLTZ_TRANSPORT':
    {
        'transport_coeff_cal'         : [int, 1, 1],             # Whether to compute transport coefficients
        'effective_mass_cal'          : [int, 1, 0],             # Whether to compute effective masses alongside transport data
        'transport_method'            : [str, 1, 'CRTA'],        # Transport approximation, e.g. constant-relaxation-time approximation
        'electron_num'                : [int, 1, None],          # Total electron count used to locate the chemical potential
        'grid'                        : [int, 3, None],          # Brillouin-zone sampling grid
        'delta_mu_range'              : [float, 2, [-5.0, 5.0]], # Chemical-potential scan window relative to the reference
        'mu_step'                     : [float, 1, 0.1],         # Chemical-potential step
        'temperature_range'           : [float, 2, [300, 300]],  # Temperature scan window in K
        'temperature_step'            : [float, 1, 50],          # Temperature step in K
        'eta'                         : [float, 1, 0.1],         # Broadening used in transport integrals
        'relax_time'                  : [float, 1, 10],          # Relaxation time in fs for CRTA-like models
        'def_pot'                     : [float, 1, 2],           # Deformation potential used by selected transport models
        'young_mod'                   : [float, 1, 240]          # Young's modulus used by selected transport models
    },

    # other utilities
    'REDUCE_BASIS':
    {
        'e_range'                     : [float, 2, None],    # Energy window used to select the reduced basis
        'threshold'                   : [float, 1, 0.02],    # Screening threshold for keeping basis components
        'band_index_range'            : [int, 2, None],      # Optional band-index window used in the reduction
        'kpoint_mode'                 : [str, 1, None]       # K-point sampling used to evaluate the reduction criterion
    },

}


# These parameters expand into a nested option dictionary.
parameter_options = {
    'package' : package,
    'kpoint_mode' : kpoint_mode,
    'integrate_mode' : integrate_mode,
    'cal_surface_method' : cal_surface_method
}

# These parameters repeat across multiple groups in the Input file.
parameter_multigroups = {
    'high_symmetry_kpoint' : 'kpoint_num',
    'kpoint_direct_coor' : 'kpoint_num',
    'lattice_vector' : None,
}

# These parameters derive shape or multiplicity from other parameters.
def operate_HR_route(nspin):
    if nspin == 1 or nspin == 4:
        return [str, 1, None]   # One HR path for spinless or spinor Hamiltonians
    elif nspin == 2:
        return [str, 2, None]   # Separate HR paths for spin-up and spin-down channels
    else:
        raise KeyError('nspin parameters setting is incorrect!')
    
def operate_kpoint_label(kpoint_num):
    return [str, kpoint_num, ['X']*kpoint_num]   # Auto-fill placeholder labels for line-mode k points

def operate_polarization_atom_type(atom_type):
    return [int, atom_type, None]                # One valence-electron count per atomic species

parameter_dependence = {
    'HR_route'                        : [['nspin'], operate_HR_route],
    'kpoint_label'                    : [['kpoint_num'], operate_kpoint_label],
    'valence_e'                       : [['atom_type'], operate_polarization_atom_type]
}
