#ifndef VELOCITY_SOLVER_H
#define VELOCITY_SOLVER_H

#include "base_data.h"
#include "xr_operation.h"
#include "band_structure_solver.h"

class velocity_solver
{
public:
    static MatrixXcd cal_velocity_1k_base(
        base_data &Base_Data,
        const VectorXcd &exp_ikR, 
        const VectorXd &eigenvalues, 
        const MatrixXcd &eigenvectors, 
        const int &alpha
    );

    static void get_velocity_matrix_alpha(
        base_data &Base_Data,
        const MatrixXd &k_direct_coor, 
        const int &alpha, 
        std::vector<MatrixXcd> &velocity_matrix
    );

    static void get_velocity_matrix(
        base_data &Base_Data,
        const MatrixXd &k_direct_coor, 
        std::array<std::vector<MatrixXcd>, 3> &velocity_matrix
    );

};

#endif