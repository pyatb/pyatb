#ifndef BERRY_CURVATURE_SOLVER_H
#define BERRY_CURVATURE_SOLVER_H

#include "base_data.h"
#include "band_structure_solver.h"
#include "velocity_solver.h"
#include "linear_response.h"

class berry_curvature_solver
{
public:
    static MatrixXd get_total_bc_fermi(
        base_data &Base_Data, 
        const MatrixXd &k_direct_coor, 
        const double &fermi_energy, 
        const int mode
    );

    static MatrixXd get_total_bc_occupiedNumber(
        base_data &Base_Data, 
        const MatrixXd &k_direct_coor, 
        const int &occupied_band_num, 
        const int mode
    );

    static void cal_total_bc_byKubo_1k(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &occupied_band_num, 
        const VectorXd &eigenvalues, 
        const MatrixXcd &eigenvectors, 
        double &total_bc_x, 
        double &total_bc_y, 
        double &total_bc_z
    );

    static void cal_total_bc_byDirect_1k(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const VectorXd &occupation_fn,
        double &total_bc_x, 
        double &total_bc_y, 
        double &total_bc_z
    );


private:
    static const Matrix3i index_m;

};


#endif