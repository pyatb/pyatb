"""
If the parameter is None, it indicates that this parameter has no default value and must be set.
If the parameter is [], it indicates that this parameter has no default value and don't need to be set in some conditions.
"""

# function switch
function_switch = {
    'INPUT_PARAMETERS'        : False,
    'LATTICE'                 : False,

    # fermi
    'BAND_STRUCTURE'          : False,
    'BANDUNFOLDING'           : False,
    'FAT_BAND'                 : False,
    'FERMI_ENERGY'            : False,
    'FERMI_SURFACE'           : False,
    'FIND_NODES'              : False,
    'JDOS'                    : False,
    'PDOS'                    : False,
    'REDUCE_BASIS'            : False,
    'SPIN_TEXTURE'            : False,
    'SURFACE_STATE'           : False,

    # berry
    'AHC'                     : False,
    'BERRY_CURVATURE'         : False,
    'BERRY_CURVATURE_DIPOLE'  : False,
    'CHERN_NUMBER'            : False,
    'CHIRALITY'               : False,
    'OPTICAL_CONDUCTIVITY'    : False,
    'POLARIZATION'            : False,
    'SHIFT_CURRENT'           : False,
    'WILSON_LOOP'             : False,
    'CPGE'                 :False,
    'DRUDE_WEIGHT'           :False
}

# these block name all have default parameters
block_can_be_empty = []

# these block function require rR matrix
need_rR_matrix = [
    'AHC',
    'BERRY_CURVATURE',
    'BERRY_CURVATURE_DIPOLE',
    'CHERN_NUMBER',
    'CHIRALITY',
    'OPTICAL_CONDUCTIVITY',
    'POLARIZATION',
    'SHIFT_CURRENT',
    'WILSON_LOOP',
    'CPGE',
    'DRUDE_WEIGHT'
]

# Some function needs this dictionary that saves k point information from the Input file
kpoint_mode = {
    'mp' : 
    {
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'mp_grid'                     : [int, 3, None]
    },

    'line': 
    {
        'kpoint_num'                  : [int, 1, None],
        'high_symmetry_kpoint'        : [float, None, 4, None],
      # 'kpoint_num_in_line'          : [int, None, None]
    },

    'direct' : 
    {
        'kpoint_num'                  : [int, 1, None],
        'kpoint_direct_coor'          : [float, None, 3, None]
    }
}

cal_surface_method = {
    'direct_diag' : 
    {
        'slab_layers'                 : [int, 1, None]
    },

    'direct_green' : 
    {
        'slab_layers'                 : [int, 1, None],
        'green_eta'                   : [float, 1, 0.001]
    },

    'green_fun' : 
    {
        'green_eta'                   : [float, 1, 0.001],
    }
}

integrate_mode = {
    'Grid' : 
    {
        'integrate_grid'              : [int, 3, [4, 4, 4]],
        'adaptive_grid'               : [int, 3, [4, 4, 4]],
        'adaptive_grid_threshold'     : [float, 1, 50.0]
    },

    'Adaptive' :
    {
        'relative_error'              : [float, 1, 1e-6],
        'absolute_error'              : [float, 1, 0.1],
        'initial_grid'                : [int, 3, [1, 1, 1]]
    }
}

INPUT = {
    'INPUT_PARAMETERS' : 
    {
        'nspin'                       : [int, 1, None],
        'package'                     : [str, 1, 'ABACUS'],
        'fermi_energy'                : [str, 1, 'Auto'],
        'fermi_energy_unit'           : [str, 1, 'eV'],
        'HR_route'                    : [str, None, None],
        'SR_route'                    : [str, 1, None],
        'rR_route'                    : [str, 1, []],
        'binary'                      : [int, 1, 0],
        'HR_unit'                     : [str, 1, 'Ry'],
        'rR_unit'                     : [str, 1, 'Bohr'],
        'max_kpoint_num'              : [int, 1, 8000],
        'sparse_format'               : [int, 1, 0]
    },

    'LATTICE' : 
    {
        'lattice_constant'            : [float, 1, None],
        'lattice_constant_unit'       : [str, 1, 'Bohr'],
        'lattice_vector'              : [float, 3, 3, None]
    },

    # fermi

    'BAND_STRUCTURE' : 
    {
        'wf_collect'                  : [int, 1, False],
        'kpoint_mode'                 : [str, 1, None]
    },

    'BANDUNFOLDING' : 
    {
        'stru_file'                   : [str, 1, None],
        'ecut'                        : [float, 1, 10],
        'band_range'                  : [int, 2, None],
        'm_matrix'                    : [float, 9, None],
        'kpoint_mode'                 : [str, 1, None]
    },

    'FAT_BAND' : 
    {
        'band_range'                  : [int, 2, None],
        'stru_file'                   : [str, 1, None],
        'kpoint_mode'                 : [str, 1, None]
    },

    'FERMI_ENERGY' : 
    {
        'temperature'                 : [float, 1, 0.0],
        'electron_num'                : [int, 1, None],
        'grid'                        : [int, 3, [10, 10, 10]],
        'epsilon'                     : [float, 1, 1e-3]
    },

    'FIND_NODES':
    {
        'energy_range'                : [float, 2, [0.0, 0.0]],
        'initial_grid'                : [int, 3, [10, 10, 10]],
        'initial_threshold'           : [float, 1, 0.1],
        'adaptive_grid'               : [int, 3, [20, 20, 20]],
        'adaptive_threshold'          : [float, 1 , 1e-3],
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]]
    },

    'BERRY_CURVATURE' : 
    {
        'method'                      : [int, 1, 0],
        'occ_band'                    : [int, 1, -1],
        'kpoint_mode'                 : [str, 1, None]
    },

    'CHERN_NUMBER' : 
    {
        'method'                      : [int, 1, 0],
        'occ_band'                    : [int, 1, -1],
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        'integrate_mode'              : [str, 1, None]
    },

    'AHC' : 
    {
        'method'                      : [int, 1, 0],
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_mode'              : [str, 1, None]
    },

    'OPTICAL_CONDUCTIVITY' : 
    {
        'occ_band'                    : [int, 1, None],
        'omega'                       : [float, 2, None],
        'domega'                      : [float, 1, None],
        'eta'                         : [float, 1, 0.01],  # unit is eV
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None],
        'method'                      : [int, 1, 1]
    },

    'POLARIZATION' : 
    {
        'occ_band'                    : [int, 1, None],
        'nk1'                         : [int, 1, 8],
        'nk2'                         : [int, 1, 8],
        'nk3'                         : [int, 1, 8],
        'atom_type'                   : [int, 1, None],
        'stru_file'                   : [str, 1, None],
        'valence_e'                   : [int, None, None]
    },

    'WILSON_LOOP' : 
    {
        'occ_band'                    : [int, 1, None],
        'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        'nk1'                         : [int, 1, 100],
        'nk2'                         : [int, 1, 100],
    },

    'SHIFT_CURRENT' : 
    {
        'occ_band'                    : [int, 1, None],
        'omega'                       : [float, 2, None],
        'domega'                      : [float, 1, None],
        'smearing_method'             : [int, 1, 1], # 0: no smearing, 1: Gauss smearing, 2: adaptive smearing
        'eta'                         : [float, 1, 0.01],  # unit is eV
        # 'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        # 'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        # 'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        # 'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'grid'                        : [int, 3, None],
        'method'                      : [int, 1, 1]
    },

    'FERMI_SURFACE' : 
    {
        'bar'                         : [float, 1, 1e-3],
        'nbands'                      : [int, 2, [0,0]],
        'kpoint_mode'                 : [str, 1, None]
    },

    'JDOS' : 
    {
        'occ_band'                    : [int, 1, None],
        'omega'                       : [float, 2, None],
        'domega'                      : [float, 1, None],
        'eta'                         : [float, 1, 0.01],  # unit is eV
        'grid'                        : [int, 3, None]
    },
    
    'PDOS' : 
    {
        'stru_file'                   : [str, 1, None],
        'e_range'                     : [float, 2, None],
        'de'                          : [float, 1, 0.01],
        'sigma'                       : [float, 1, 1e-3],
        'kpoint_mode'                 : [str, 1, None]
    },

    'SPIN_TEXTURE' : 
    {
        'nband'                       : [int, 1, None],
        'kpoint_mode'                 : [str, 1, None]
    },

    'BERRY_CURVATURE_DIPOLE':
    {
        'omega'                       : [float, 2, None],
        'domega'                      : [float, 1, None],
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_mode'              : [str, 1, None]
    },

    'CHIRALITY':
    {
        'method'                      : [int, 1, 0],
        'k_vect'                      : [float, 3, [0.0, 0.0, 0.0]],
        'radius'                      : [float, 1, 0.01],
        'point_num'                   : [int, 1, 1000]
    },

    'SURFACE_STATE':
    {
        'cal_surface_method'          : [str, 1, None],
        'surface_direction'           : [str, 1, 'c'],
        'energy_windows'              : [float, 3, None],
        'coupling_layers'             : [int, 1, None],
        'kpoint_mode'                 : [str, 1, None]
    },

    'CPGE':
    {
        'omega'                       : [float, 2, None],
        'domega'                      : [float, 1, None],
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_mode'              : [str, 1, None]
    },

    'DRUDE_WEIGHT':
    {
        'omega'                       : [float, 2, None],
        'domega'                      : [float, 1, None],
        #'k_start'                     : [float, 3, [0.0, 0.0, 0.0]],
        #'k_vect1'                     : [float, 3, [1.0, 0.0, 0.0]],
        #'k_vect2'                     : [float, 3, [0.0, 1.0, 0.0]],
        #'k_vect3'                     : [float, 3, [0.0, 0.0, 1.0]],
        'integrate_mode'              : [str, 1, None]
    },

    'REDUCE_BASIS':
    {
        'e_range'                     : [float, 2, None],
        'threshold'                   : [float, 1, 0.02],
        'band_index_range'            : [int, 2, None],
        'kpoint_mode'                 : [str, 1, None]
    },

}


# these parameters have many options
parameter_options = {
    'kpoint_mode' : kpoint_mode,
    'integrate_mode' : integrate_mode,
    'cal_surface_method' : cal_surface_method
}

# these parameters are displayed in multiple groups in the Input file
parameter_multigroups = {
    'high_symmetry_kpoint' : 'kpoint_num',
    'kpoint_direct_coor' : 'kpoint_num',
    'lattice_vector' : None,
}

# these parameters depend on other parameter
def operate_HR_route(nspin):
    if nspin == 1 or nspin == 4:
        return [str, 1, None]
    elif nspin == 2:
        return [str, 2, None]
    else:
        raise KeyError('nspin parameters setting is incorrect!')

def operate_polarization_atom_type(atom_type):
    return [int, atom_type, None]

parameter_dependence = {
    'HR_route'                        : [['nspin'], operate_HR_route],
    'valence_e'                       : [['atom_type'], operate_polarization_atom_type]
}
