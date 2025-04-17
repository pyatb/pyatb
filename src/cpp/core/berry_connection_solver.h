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

    static MatrixXcd get_berry_connection_sumOver(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &alpha
    ); 

    static void get_berry_connection_sumOver_alldirection(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        std::array<MatrixXcd, 3> &berry_connection
    ); 

    static MatrixXcd get_partial_berry_connection(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &alpha,
        const int &beta
    );

    static void get_rnm_and_drnm(
        base_data &Base_Data,
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        std::array<MatrixXcd, 3> &r_nm,
        std::array<MatrixXcd, 9> &d_r_nm,
        const double &degenerate_eta = 0.04
    ); 

};

#endif