import numpy as np
import re
from pyatb.constants import Ang_to_Bohr

class Atom:
    def __init__(self, species):
        self.species = species

    def set_positions_car(self, positions_car:np.ndarray):
        self.atom_num = positions_car.shape[0]
        self.cartesian_coor = positions_car

    def set_numerical_orb(self, orb_file):
        self.orb_file = orb_file

    def read_numerical_orb(self):
        with open(self.orb_file, 'r') as f:
            line = f.readline()
            line = f.readline()
            line = f.readline()
            line = f.readline()
            self.l_max = int(f.readline().split()[-1])
            self.orbital_num = list()
            for i in range(self.l_max+1):
                line = f.readline().split()
                self.orbital_num.append(int(line[-1]))
            line = f.readline()
            line = f.readline()
            line = f.readline()
            self.mesh = int(f.readline().split()[-1])
            self.dr = float(f.readline().split()[-1])
            total_orb_num = np.sum(np.array(self.orbital_num))
            self.orbit = np.zeros([total_orb_num, self.mesh], dtype=float)
            line_num, remain_num = divmod(self.mesh, 4)
            for i in range(total_orb_num):
                line = f.readline()
                line = f.readline()
                count = 0
                for j in range(line_num):
                    line = f.readline().split()
                    for k in range(4):
                        self.orbit[i, count] = float(line[k])# / Ang_to_Bohr
                        count += 1
                if remain_num:
                    line = f.readline().split()
                    for j in range(remain_num):
                        self.orbit[i, count] = float(line[j])# / Ang_to_Bohr
                        count += 1


def skip_notes(line):
    """ remove the comment on a line in the file. """
    line = re.compile(r"#.*").sub("", line)
    line = re.compile(r"//.*").sub("", line)
    line = line.strip()
    return line

def read_stru(stru_filename='STRU'):
    with open(stru_filename, 'r') as f:
        raw_file = list()
        for line in f:
            line = skip_notes(line)
            if line:
                raw_file.append(line)

    ATOMIC_SPECIES = 'ATOMIC_SPECIES'
    NUMERICAL_ORBITAL = 'NUMERICAL_ORBITAL'
    LATTICE_CONSTANT = 'LATTICE_CONSTANT'
    LATTICE_VECTORS = 'LATTICE_VECTORS'
    ATOMIC_POSITIONS = 'ATOMIC_POSITIONS'

    try:
        ATOMIC_SPECIES_index = raw_file.index(ATOMIC_SPECIES)
        NUMERICAL_ORBITAL_index = raw_file.index(NUMERICAL_ORBITAL)
        LATTICE_CONSTANT_index = raw_file.index(LATTICE_CONSTANT)
        LATTICE_VECTORS_index = raw_file.index(LATTICE_VECTORS)
        ATOMIC_POSITIONS_index = raw_file.index(ATOMIC_POSITIONS)
    except:
        print('The structure file format is incorrect !')

    atom_type_num = NUMERICAL_ORBITAL_index - ATOMIC_SPECIES_index - 1
    stru_atom = [None] * atom_type_num
    for i in range(atom_type_num):
        temp_species = raw_file[ATOMIC_SPECIES_index+i+1].split()[0]
        stru_atom[i] = Atom(temp_species)
        temp_orb = raw_file[NUMERICAL_ORBITAL_index+i+1]
        stru_atom[i].set_numerical_orb(temp_orb)

    A = float(raw_file[LATTICE_CONSTANT_index+1]) / Ang_to_Bohr
    V = np.zeros((3, 3), dtype=float)
    for i in range(3):
        temp_v = raw_file[LATTICE_VECTORS_index+i+1].split()
        for j in range(3):
            V[i, j] = float(temp_v[j])
    
    if raw_file[ATOMIC_POSITIONS_index+1] == 'Direct':
        is_car = False
    else:
        is_car = True

    count = 2
    for i in range(atom_type_num):
        atom_type = raw_file[ATOMIC_POSITIONS_index+count]
        if atom_type != stru_atom[i].species:
            raise ValueError('The structure file format is incorrect !')
        count += 2
        atom_num = int(raw_file[ATOMIC_POSITIONS_index+count])
        atom_postion_car = np.zeros([atom_num, 3], dtype=float)
        count += 1
        for ia in range(atom_num):
            temp_positon = raw_file[ATOMIC_POSITIONS_index+count+ia].split()
            atom_postion_car[ia, 0] = float(temp_positon[0])
            atom_postion_car[ia, 1] = float(temp_positon[1])
            atom_postion_car[ia, 2] = float(temp_positon[2])
        count += atom_num
        if is_car:
            stru_atom[i].set_positions_car(atom_postion_car * A)
        else:
            stru_atom[i].set_positions_car(atom_postion_car @ V * A)

    # print('atom_type_num = ', atom_type_num)

    # for i in stru_atom:
    #     i.read_numerical_orb()

    return stru_atom
