#ifndef BERRY_PHASE_SOLVER_H
#define BERRY_PHASE_SOLVER_H

#include "base_data.h"
#include "xr_operation.h"
#include "band_structure_solver.h"

class berry_phase_solver
{
public:
    // k point is direct coordinate
    static double get_berry_phase_of_loop(base_data &Base_Data, const MatrixXd &k_direct_coor_loop, const int &occupied_band_num);

    static VectorXd calculate_wilson_loop(base_data &Base_Data, const MatrixXd &k_direct_coor_loop, const int &occupied_band_num);
    
};


#endif