#ifndef SHIFT_CURRENT_SOLVER_H
#define SHIFT_CURRENT_SOLVER_H

#include "base_data.h"
#include "band_structure_solver.h"
#include "berry_connection_solver.h"
#include "velocity_solver.h"
#include "linear_response.h"
#include <set>
#include <map>

class shift_current_solver
{
public:

    void set_parameters(
        const int &nspin,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const int &smearing_method,
        const double &eta
    );

    MatrixXd get_shift_current_conductivity(
        base_data &Base_Data,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num,
        const int &occupied_num,
        const int &method
    );

private:

    MatrixXd get_shift_current_conductivity_ik(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const int &occupied_num,
        const int &method
    );

    // m is the band indicator of the occupied state, 
    // n is the band indicator of the unoccupied state, 
    // and n, m is a pair of indicators, so the dimensions of m and n are the same.
    MatrixXd get_Inm_ik_direct(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const std::vector<int> m,
        const std::vector<int> n
    );

    // m is the band indicator of the occupied state, 
    // n is the band indicator of the unoccupied state, 
    // and n, m is a pair of indicators, so the dimensions of m and n are the same.
    MatrixXd get_Inm_ik_sumOver(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const std::vector<int> m,
        const std::vector<int> n
    );

    int nspin = 1;
    int omega_num;
    double domega;
    double start_omega;
    int smearing_method = 1;
    double eta;
    
};

#endif