#ifndef POCKELS_SOLVER_H
#define POCKELS_SOLVER_H

#include "base_data.h"
#include "band_structure_solver.h"
#include "berry_connection_solver.h"
#include "velocity_solver.h"
#include "linear_response.h"
#include <set>
#include <map>

class pockels_solver
{
public:

    void set_parameters(
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &fermi_energy,
        const double &omega1
    );

    MatrixXcd get_pockels(
        base_data &Base_Data,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num
    );

private:

    MatrixXcd get_pockels_ik(
        base_data &Base_Data,
        const VectorXcd &exp_ikR
    );

    // m is the band indicator of the occupied state, 
    // n is the band indicator of the unoccupied state, 
    // and n, m is a pair of indicators, so the dimensions of m and n are the same.
    

    int omega_num;
    double domega;
    double start_omega;
    double fermi_energy;
    double omega1;
    
};

#endif