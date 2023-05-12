#ifndef BERRY_CONNECTION_SOLVER_H
#define BERRY_CONNECTION_SOLVER_H

#include "base_data.h"
#include "xr_operation.h"
#include "linear_response.h"

class berry_connection_solver
{
public:
    static MatrixXcd get_berry_connection(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &alpha
    );

    static MatrixXcd get_partial_berry_connection(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &alpha,
        const int &beta
    );
};

#endif