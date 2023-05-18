#ifndef BAND_STRUCTURE_SOLVER_H
#define BAND_STRUCTURE_SOLVER_H


#include "base_data.h"
#include "xr_operation.h"
#include "tools.h"

class band_structure_solver
{
public:
    static void get_eigenvalues_1k(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        VectorXd &eigenvalues
    );
    
    static void get_eigenvalues_eigenvectors_1k(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        VectorXd &eigenvalues, 
        MatrixXcd &eigenvectors
    );

    static void get_eigenvalues(
        base_data &Base_Data, 
        const MatrixXd &k_direct_coor, 
        std::vector<VectorXd> &eigenvalues
    );
    
    static void get_eigenvalues_eigenvectors(
        base_data &Base_Data, 
        const MatrixXd &k_direct_coor, 
        std::vector<VectorXd> &eigenvalues, 
        std::vector<MatrixXcd> &eigenvectors
    );
};

#endif