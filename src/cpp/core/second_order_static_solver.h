#ifndef SECOND_ORDER_STATIC_SOLVER_H
#define SECOND_ORDER_STATIC_SOLVER_H

#include "base_data.h"
#include "band_structure_solver.h"
#include "berry_connection_solver.h"
#include "velocity_solver.h"
#include "linear_response.h"
#include <set>
#include <map>

class second_order_static_solver
{
public:
    MatrixXcd get_second_order_static(
        base_data &Base_Data,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num,
        const float &fermi_energy
    );

private:

    MatrixXcd get_static_ik(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const float &fermi_energy
    );

    // m is the band indicator of the occupied state, 
    // n is the band indicator of the unoccupied state, 
    // and n, m is a pair of indicators, so the dimensions of m and n are the same.
    
    
};

#endif