from pyatb import RANK, OUTPUT_PATH, timer
from pyatb.io.default_input import function_switch, block_can_be_empty, need_rR_matrix
from pyatb.io.default_input import INPUT, parameter_options, parameter_multigroups, parameter_dependence
from pyatb.constants import Ry_to_eV, Ang_to_Bohr

import re
import json
from collections import defaultdict
import numpy as np
import copy



def skip_notes(line: str):
    """ remove the comment on a line in the file. """
    line = re.compile(r"#.*").sub("", line)
    line = re.compile(r"//.*").sub("", line)
    line = line.strip()
    return line

def get_file_block(input_filename: str) -> dict:
    """
    Get the All block of the input_filename file.

    Explanation:
        The block refers to the following form in the input_filename file:
            '''
            INPUT_PARAMETERS
            {
                nspin          1
                package        ABACUS
            }
            '''
        The INPUT_PARAMETERS is block name, the content is the parameter list in the block INPUT_PARAMETERS{}.

    Return:
        The return value file_block is a dict, the form is {'INPUT_PARAMETERS' : ['nspin', '1', 'package', 'ABACUS'], etc.}
    """
    with open(input_filename, 'r') as f:
        raw_file = ''
        for line in f:
            line = skip_notes(line)
            raw_file += line + '\n'
        
        processed_file = re.compile(r'[_\w]+\s*\{[\d\D]*?\}').findall(raw_file)
        file_block_tem = [None for i in range(len(processed_file))]
        compile = re.compile(r'[\{\},\s]')
        for i in range(len(processed_file)):
            file_block_tem[i] = compile.split(processed_file[i])
            file_block_tem[i] = list(filter(None, file_block_tem[i]))
    
    file_block = {i_block[0] : i_block[1:] for i_block in file_block_tem}

    return file_block

def get_block_data(block_name: str, file_block: dict): 
    """ Get the data in the block_name from file_block dict. """
    can_be_empty = False  # determine whether all block parameters have default values, True -> Yes, False -> No
    if block_name in block_can_be_empty:
        can_be_empty = True

    if block_name in file_block:
        data = file_block.get(block_name)
        if data:
            function_switch[block_name] = True
            return data
        elif can_be_empty:
            function_switch[block_name] = True
            return None
        else:
            raise ValueError(block_name + ' is empty!')

def get_general_parameter(parameter_name: str, default: None, data: None):
    """ 
    Extract information with the parameter_name from data.

    Parameters:
        parameter_name : parameter name, such as 'nspin'.
        default : describe the type, size, and default value of the parameter, such as [str, 1, 'Auto'], or [str, 1, None].
        data : save the list of parametes, such as ['nspin', '1', 'grid', '4', '4', '4'].
    """
    for index, value in enumerate(data):
        if value == parameter_name:
            if default[1] == 1:
                return default[0](data[index+1])
            else:
                tem = []
                for i_size in range(default[1]):
                    tem.append(default[0](data[index + 1 + i_size]))
                return tem

    if default[2] is None:
        raise KeyError('key parameter missing: ' + parameter_name)
    else:
        return default[2]

def get_multiline_parameters(parameter_name, default, data):
    """ 
    Extract information with the parameter_name from data.

    Parameters:
        parameter_name : parameter name, such as 'lattice_vector'.
        default : describe the type, group number, size in every group, and default value of the parameter, such as [float, 3， 3, None].
        data : save the list of parametes, such as ['lattice_vector', '1', '0', '0', '0', '1', '0', '0', '0', '1'].
    """
    for index, value in enumerate(data):
        if value == parameter_name:
            tem = []
            count = index + 1
            for i_group in range(default[1]):
                group_list = []
                for j_size in range(default[2]):
                    group_list.append(default[0](data[count]))
                    count = count + 1
                tem.append(group_list)
            return tem

    if default[3] is None:
        raise KeyError('key parameter missing: ' + parameter_name)
    else:
        return default[3]

def update_INPUT(input_filename):
    """ Update the dictionary INPUT in default_input.py according to the input_filename file. """
    file_block = get_file_block(input_filename)

    for block_name in INPUT.keys():
        block_data = get_block_data(block_name, file_block)
        if block_data:
            block_parameters = INPUT[block_name]

            # search optional parameters
            for i in parameter_options.keys():
                if i in block_parameters.keys():
                    block_parameters[i][-1] = get_general_parameter(i, block_parameters[i], block_data)
                    option = block_parameters[i][-1]
                    option_dict = copy.deepcopy(parameter_options[i][option])
                    block_parameters.update(option_dict)

            # search general parameters
            for i in block_parameters.keys():
                if i not in parameter_multigroups and i not in parameter_dependence:
                    block_parameters[i][-1] = get_general_parameter(i, block_parameters[i], block_data)

            # search multigroups parameters
            for i in parameter_multigroups.keys():
                if i in block_parameters.keys():
                    if block_parameters[i][1] is None:
                        block_parameters[i][1] = block_parameters[parameter_multigroups[i]][-1]
                    block_parameters[i][-1] = get_multiline_parameters(i, block_parameters[i], block_data)

            # search dependent parameters
            for i in parameter_dependence.keys():
                if i in block_parameters.keys():
                    tem_list = []
                    for tem_value in parameter_dependence[i][0]:
                        tem_list.append(block_parameters[tem_value][-1])
                    tem_fun = parameter_dependence[i][1]
                    block_parameters[i] = tem_fun(*tem_list)
                    block_parameters[i][-1] = get_general_parameter(i, block_parameters[i], block_data)
                    

    # rearrange the parameters in the dictionary INPUT
    for block_name in INPUT.keys():
        block_parameters = INPUT[block_name]
        for i in block_parameters.keys():
            block_parameters[i] = block_parameters[i][-1]

def parameter_require_additional_operations():
    # for fermi_energy
    tem_value = INPUT['INPUT_PARAMETERS']['fermi_energy']
    if tem_value != 'Auto':
        tt_value = INPUT['INPUT_PARAMETERS']['fermi_energy_unit']
        if tt_value == 'Ry':
            INPUT['INPUT_PARAMETERS']['fermi_energy'] = float(tem_value) * Ry_to_eV
            INPUT['INPUT_PARAMETERS']['fermi_energy_unit'] = 'eV'
        elif tt_value == 'eV':
            INPUT['INPUT_PARAMETERS']['fermi_energy'] = float(tem_value)

    # for high_symmetry_kpoint
    for block_name in INPUT.keys():
        block_parameters = INPUT[block_name]
        high_symmetry_kpoint = block_parameters.get('high_symmetry_kpoint')
        if high_symmetry_kpoint:
            block_parameters['kpoint_num_in_line'] = [int(i.pop()) for i in high_symmetry_kpoint]

    # for lattice constant
    tem_value = INPUT['LATTICE']['lattice_constant_unit']
    if tem_value == 'Bohr':
        INPUT['LATTICE']['lattice_constant'] /= Ang_to_Bohr
        tem_value = INPUT['LATTICE']['lattice_constant_unit'] = 'Angstrom'
    elif tem_value == 'Angstrom':
        pass

def check():
    if not function_switch['INPUT_PARAMETERS']:
        raise KeyError('you have to set up INPUT_PARAMETERS block')
    
    if not function_switch['LATTICE']:
        raise KeyError('you have to set up LATTICE block')

    if INPUT['INPUT_PARAMETERS']['fermi_energy'] == 'Auto' and not function_switch['FERMI_ENERGY']:
        raise KeyError('if fermi_energy is Auto, please set FERMI_ENERGY block')

    bool_need_rR = False
    for key in need_rR_matrix:
        bool_need_rR = bool_need_rR or function_switch[key]

    if INPUT['INPUT_PARAMETERS']['package'] == 'ABACUS' and INPUT['INPUT_PARAMETERS']['rR_route'] == []:
        if bool_need_rR:
            raise KeyError('you have to set up rR_route parameters in INPUT_PARAMETERS block')

def print_json(dict, filename):
    import os
    with open(os.path.join(OUTPUT_PATH, filename), 'w') as f:
        json.dump(obj=dict, fp=f, indent=4)

def read_input(input_filename):
    # start time
    timer.start('read_input', 'Read Input File')

    update_INPUT(input_filename)

    parameter_require_additional_operations()

    check()

    if RANK == 0:
        print_json(INPUT, "input.json")

    new_INPUT = dict()
    for function, switch in function_switch.items():
        if switch:
            new_INPUT[function] = INPUT[function]

    # list to numpy.array
    for block_name in new_INPUT.keys():
        block_parameters = new_INPUT[block_name]
        for key, values in block_parameters.items():
            if isinstance(values, list):
                block_parameters[key] = np.array(values)

    bool_need_rR = False
    for key in need_rR_matrix:
        bool_need_rR = bool_need_rR or function_switch[key]

    # end time
    timer.end('read_input', 'Read Input File')

    return new_INPUT, function_switch, bool_need_rR
    
    

    
