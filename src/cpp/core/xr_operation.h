#ifndef XR_OPERATION
#define XR_OPERATION

#include "base_data.h"


class xr_operation
{
public:
    static MatrixXcd get_Hk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR
    );

    static MatrixXcd get_Sk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR
    );

    static MatrixXcd get_partial_Hk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &partial_alpha
    );

    static MatrixXcd get_partial_Sk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &partial_alpha
    );

    static MatrixXcd get_partial2_Hk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &partial_alpha,
        const int &partial_beta
    );

    static MatrixXcd get_partial2_Sk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &partial_alpha,
        const int &partial_beta
    );

    static MatrixXcd get_rk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &alpha
    );

    static MatrixXcd get_rk_S(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const MatrixXcd &partial_alpha_Sk, 
        const int &alpha
    );

    static MatrixXcd get_partial_rk(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &partial_alpha, 
        const int &beta
    );

    static MatrixXcd get_partial_rk_S(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const MatrixXcd &partial_alpha_beta_Sk, 
        const int &partial_alpha, 
        const int &beta
    );

    static MatrixXcd get_Omega_extra(
        base_data &Base_Data, 
        const VectorXcd &exp_ikR, 
        const int &alpha, 
        const int &beta
    );

    static MatrixXcd inner_product_twoPoints(
        base_data &Base_Data, 
        const VectorXd &k_direct_coor_start, 
        const VectorXd &k_direct_coor_end
    );

    static std::vector<MatrixXcd> inner_product_along_loop(
        base_data &Base_Data, 
        const MatrixXd &k_direct_coor_loop
    );

    static void Fourier_R_to_k_SurfaceState_bulk_Xk(
        base_data &Base_Data,
        const int &direction,
        const MatrixXd &k_direct_coor,
        const int &coupling_layers,
        std::vector<MatrixXcd> &bulk_Xk,
        const MatrixXcd &XR_upperTriangleOfDenseMatrix
    );

    static void SurfaceState_Xk(
        base_data &Base_Data,
        const int &direction,
        const MatrixXd &k_direct_coor,
        const int &coupling_layers,
        const char &X,
        std::vector<MatrixXcd> &Xk00,
        std::vector<MatrixXcd> &Xk01
    );

};


#endif