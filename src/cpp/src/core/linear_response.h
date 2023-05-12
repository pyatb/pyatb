#ifndef LINEAR_RESPONSE_H
#define LINEAR_RESPONSE_H

#include "base_data.h"
#include "xr_operation.h"
#include "band_structure_solver.h"
#include "tools.h"

class linear_response
{
public:
    // static double get_partial_eigenvalue(
    //     const MatrixXcd &H_alpha,
    //     const MatrixXcd &S_alpha,
    //     const double &eigenvalue,
    //     const VectorXcd &eigenvector
    // );

    static double get_partial_eigenvalue(
        const MatrixXcd &H_alpha,
        const MatrixXcd &S_alpha,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &band_index
    );

    static double get_partial2_eigenvalue(
        const MatrixXcd &H_alpha,
        const MatrixXcd &H_beta_alpha,
        const MatrixXcd &S_alpha,
        const MatrixXcd &S_beta_alpha,
        const double &eigenvalue,
        const double &eigenvalue_beta,
        const VectorXcd &eigenvector,
        const VectorXcd &eigenvector_beta
    );


    static VectorXcd get_partial_eigenvectors_1k_Sternheimer(
        const MatrixXcd &H,
        const MatrixXcd &S,
        const MatrixXcd &H_alpha,
        const MatrixXcd &S_alpha,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &band_index
    );

    static VectorXcd get_partial_eigenvectors_1k_SumOverStates(
        const MatrixXcd &H_alpha,
        const MatrixXcd &S_alpha,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &band_index
    );


    static MatrixXcd get_partial_eigenvectors_1k(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &alpha,
        const int &mode
    );

    static VectorXcd get_partial2_eigenvectors_1k_Sternheimer(
        const MatrixXcd &H,
        const MatrixXcd &S,
        const MatrixXcd &H_alpha,
        const MatrixXcd &S_alpha,
        const MatrixXcd &H_beta,
        const MatrixXcd &S_beta,
        const MatrixXcd &H_beta_alpha,
        const MatrixXcd &S_beta_alpha,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const VectorXcd &eigenvector_alpha,
        const VectorXcd &eigenvector_beta,
        const int &band_index
    );

    static VectorXcd get_partial2_eigenvectors_1k_SumOverStates(
        const MatrixXcd &H_alpha,
        const MatrixXcd &H_beta,
        const MatrixXcd &H_beta_alpha,
        const MatrixXcd &S,
        const MatrixXcd &S_alpha,
        const MatrixXcd &S_beta,
        const MatrixXcd &S_beta_alpha,
        const VectorXd &eigenvalues,
        const double &E_alpha,
        const double &E_beta,
        const MatrixXcd &eigenvectors,
        const VectorXcd &eigenvector_alpha,
        const VectorXcd &eigenvector_beta,
        const int &band_index
    );

    static MatrixXcd get_partial2_eigenvectors_1k(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &alpha,
        const int &beta,
        const int &mode
    );

    static MatrixXcd get_D(
        const MatrixXcd &H_alpha_bar,
        const MatrixXcd &S_alpha_bar,
        const MatrixXcd &A_alpha_bar_dagger,
        const VectorXd &eigenvalues
    );

    static MatrixXcd get_D_degenerate(
        const MatrixXcd &H_alpha_bar,
        const MatrixXcd &S_alpha_bar,
        const MatrixXcd &A_alpha_bar_dagger,
        const VectorXd &eigenvalues,
        const double &degenerate_eta=0.04
    );

    static MatrixXcd get_D_e(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR,
        const MatrixXcd &Hk,
        const MatrixXcd &Sk,
        const VectorXd &eigenvalues,
        const MatrixXcd &eigenvectors,
        const int &alpha
    );

private:
    static const double degenerate_threshold;

};

#endif