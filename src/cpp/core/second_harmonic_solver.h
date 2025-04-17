#ifndef SECOND_HARMONIC_SOLVER_H
#define SECOND_HARMONIC_SOLVER_H

#include "base_data.h"
#include "band_structure_solver.h"
#include "berry_connection_solver.h"
#include "velocity_solver.h"
#include "linear_response.h"
#include <set>
#include <map>

class second_harmonic_solver
{
public:

    void set_parameters(
        const int &method,
        const double &eta,
        const int &omega_num,
        const double &domega,
        const double &start_omega
    );

    MatrixXcd get_second_harmonic(
        base_data &Base_Data,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num,
        const float &fermi_energy
    );

private:

    MatrixXcd get_second_harmonic_ik(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const float &fermi_energy
    );

    // m is the band indicator of the occupied state, 
    // n is the band indicator of the unoccupied state, 
    // and n, m is a pair of indicators, so the dimensions of m and n are the same.
    
    int method;
    double eta;
    int omega_num;
    double domega;
    double start_omega;
    
};

#endif