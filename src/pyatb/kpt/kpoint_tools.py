import numpy as np
import os

def print_kpoint_line_style_plot_assistance_file(
    output_path,
    high_symmetry_kpoint,
    kpoint_num_in_line,
    kpoint_label,
    lattice_constant,
    lattice_vector
):
    reciprocal_vector = np.linalg.inv(lattice_vector).transpose() * 2 * np.pi / lattice_constant
    def direct_to_cartesian(k_direct_coor):
        kvect_cartesian_coor = k_direct_coor @ reciprocal_vector
        return kvect_cartesian_coor

    with open(os.path.join(output_path, 'high_symmetry_kpoint.dat'), 'w') as f:
        x_coor = np.zeros(high_symmetry_kpoint.shape[0], dtype=float)
        for i in range(high_symmetry_kpoint.shape[0]-1):
            if  kpoint_num_in_line[i] <= 1:
                distance = 0
                kpoint_label[i] = kpoint_label[i] + '|' + kpoint_label[i+1]
                kpoint_label[i+1] = kpoint_label[i]
            else:
                distance = np.linalg.norm((direct_to_cartesian(high_symmetry_kpoint[i+1, :]) - direct_to_cartesian(high_symmetry_kpoint[i, :])))
            x_coor[i+1] = x_coor[i] + distance

        f.write('# K-label / x coordinate of band-figure / direct coordinate in the Brillouin zone / kpoint number in line\n')
        for i in range(high_symmetry_kpoint.shape[0]):
            f.write('%-5s %10.6f %10.6f %10.6f %10.6f %10d\n'%(kpoint_label[i], x_coor[i], high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))

    x_coor_array = np.array([], dtype=float)
    for i in range(high_symmetry_kpoint.shape[0]-1):
        x_coor_array = np.concatenate((x_coor_array, np.linspace(x_coor[i], x_coor[i+1], kpoint_num_in_line[i], endpoint=False)))
    x_coor_array = np.append(x_coor_array, x_coor[-1])
    np.savetxt(os.path.join(output_path, 'x_coor_array.dat'), x_coor_array, fmt='%10.6f')