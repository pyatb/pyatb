// A two-level, R=0 model: swapping lattice basis rows changes handedness,
// while Cartesian position matrix elements and the physical model stay fixed.
#include "optical_conductivity_solver.h"
#include <iomanip>

class interface_python {
public:
    static void populate(base_data &data, bool left) {
        data.lattice_constant = 1.0;
        data.lattice_vector = Matrix3d::Zero();
        data.lattice_vector.diagonal() << 2.0, 3.0, 4.0;
        if (left) data.lattice_vector.row(1).swap(data.lattice_vector.row(2));
        data.reciprocal_vector = data.lattice_vector.inverse().transpose();
        data.basis_num = 2;
        data.R_num = 1;
        data.R_direct_coor = MatrixXd::Zero(1, 3);
        data.R_cartesian_coor = MatrixXd::Zero(1, 3);
        data.HR_upperTriangleOfDenseMatrix = MatrixXcd::Zero(1, 3);
        data.HR_upperTriangleOfDenseMatrix << -1.0, 0.0, 1.0;
        data.SR_upperTriangleOfDenseMatrix = MatrixXcd::Zero(1, 3);
        data.SR_upperTriangleOfDenseMatrix << 1.0, 0.0, 1.0;
        for (int axis = 0; axis < 3; ++axis) {
            data.rR_upperTriangleOfDenseMatrix[axis] = MatrixXcd::Zero(1, 3);
            data.rR_upperTriangleOfDenseMatrix[axis](0, 1) = 1.0 / (axis + 1);
        }
    }
};

int main() {
    std::cout << std::setprecision(17);
    for (int left = 0; left < 2; ++left) {
        base_data data;
        interface_python::populate(data, left);
        for (int spin : {1, 4}) {
            for (int method : {0, 1}) {
                optical_conductivity_solver solver;
                MatrixXd kpoints = MatrixXd::Zero(1, 3);
                solver.set_parameters(spin, 61, 0.1, 0.0, 0.1, 1, kpoints, 1);
                MatrixXcd sigma = MatrixXcd::Zero(9, 61);
                MatrixXcd chi = MatrixXcd::Zero(9, 61);
                solver.get_optical_conductivity_by_kubo(data, method, sigma, chi);
                Matrix<double, 9, 1> static_chi = Matrix<double, 9, 1>::Zero();
                solver.get_static_dielectric_function_by_kubo(data, static_chi);
                for (int w = 0; w < 61; ++w) {
                    std::cout << left << " " << spin << " " << method << " " << w
                              << " " << data.get_primitiveCell_volume();
                    for (int axis : {0, 4, 8})
                        std::cout << " " << sigma(axis, w).real()
                                  << " " << chi(axis, w).real()
                                  << " " << chi(axis, w).imag()
                                  << " " << static_chi(axis);
                    std::cout << "\n";
                }
            }
        }
    }
}
