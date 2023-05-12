#ifndef BANDUNFOLDING_SOLVER_H
#define BANDUNFOLDING_SOLVER_H

#include "base_data.h"
#include "cell_atom.h"
#include "band_structure_solver.h"

class bandunfolding_solver
{
public:
    bandunfolding_solver();
    ~bandunfolding_solver();

    void set_M_matrix(
        const double &lattice_constant,
        const Matrix3d &lattice_vector,
        const Matrix3d &M_matrix
    );

    void output_spectral_weight(
        base_data &Base_Data, 
        const MatrixXd &kvect_direct, 
        const double &ecut,
        const int &min_bandindex,
        const int &max_bandindex,
        const int &nspin,
        MatrixXd &P,
        MatrixXd &E
    );

private:
    MatrixXd generate_GVectors_pw(const double &ecut);

    void generate_supercell_kpoint(
        const MatrixXd &unitcell_kvect_direct, 
        std::vector<int> &unit2super_tag,              // store the correspondence between the supercell k point and the primitive cell k point.
        std::vector<std::vector<int>> &super2unit_tag, // store the correspondence between the supercell k point and the primitive cell k point.
        MatrixXd &supercell_kvec_d
    );


    double lattice_constant;
    Matrix3d lattice_vector;
    Matrix3d reciprocal_vector;
    Matrix3d M_matrix;
    Matrix3d unitcell_a;
    Matrix3d unitcell_g;
};

#endif