#include "berry_phase_solver.h"
#include "tools.h"

double berry_phase_solver::get_berry_phase_of_loop(base_data &Base_Data, const MatrixXd &k_direct_coor_loop, const int &occupied_band_num)
{
    int kpoint_num = k_direct_coor_loop.rows();    
    complex<double> det = 1.0;
    
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor_loop);
    std::vector<MatrixXcd> eigenvectors(3);
    
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd inner_product;

        if (ik == 0)
        {
            VectorXd eigenvalues;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors[0]);
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik+1), eigenvalues, eigenvectors[1]);
            inner_product = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(ik+1));
            inner_product = eigenvectors[0].adjoint() * inner_product * eigenvectors[1];
        }
        else if (ik == kpoint_num - 1)
        {
            inner_product = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(0));
            inner_product = eigenvectors[1].adjoint() * inner_product * eigenvectors[0];
        }
        else
        {
            VectorXd eigenvalues;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik+1), eigenvalues, eigenvectors[2]);
            inner_product = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(ik+1));
            inner_product = eigenvectors[1].adjoint() * inner_product * eigenvectors[2];
            eigenvectors[1].swap(eigenvectors[2]);
        }

        complex<double> tem_det = inner_product.block(0, 0, occupied_band_num, occupied_band_num).determinant();
        det *= tem_det;
    }

    double phase = log(det).imag();
    return phase;
}

VectorXd berry_phase_solver::calculate_wilson_loop(base_data &Base_Data, const MatrixXd &k_direct_coor_loop, const int &occupied_band_num)
{
    int kpoint_num = k_direct_coor_loop.rows();
    MatrixXcd Lambda = MatrixXcd::Identity(occupied_band_num, occupied_band_num);

    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor_loop);
    std::vector<MatrixXcd> eigenvectors(3);

    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd inner_product;

        if (ik == 0)
        {
            VectorXd eigenvalues;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors[0]);
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik+1), eigenvalues, eigenvectors[1]);
            inner_product = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(ik+1));
            inner_product = eigenvectors[0].adjoint() * inner_product * eigenvectors[1];
        }
        else if (ik == kpoint_num - 1)
        {
            inner_product = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(0));
            inner_product = eigenvectors[1].adjoint() * inner_product * eigenvectors[0];
        }
        else
        {
            VectorXd eigenvalues;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik+1), eigenvalues, eigenvectors[2]);
            inner_product = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(ik+1));
            inner_product = eigenvectors[1].adjoint() * inner_product * eigenvectors[2];
            eigenvectors[1].swap(eigenvectors[2]);
        }

        JacobiSVD<MatrixXcd> svd;
        svd.compute(inner_product.block(0, 0, occupied_band_num, occupied_band_num), ComputeThinU | ComputeThinV);
        MatrixXcd U = svd.matrixU();
        MatrixXcd V = svd.matrixV();
        Lambda = V * U.adjoint() * Lambda;
    }

    ComplexEigenSolver<MatrixXcd> eigenSolver; 
    eigenSolver.compute(Lambda);
    VectorXcd Lambda_eigenvalues = eigenSolver.eigenvalues();

    VectorXd wilson_phase;
    wilson_phase.setZero(occupied_band_num);
    for (int ib = 0; ib < occupied_band_num; ++ib)
    {
        wilson_phase(ib) = log(Lambda_eigenvalues[ib]).imag() / TWO_PI;
        wilson_phase(ib) = wilson_phase[ib] - floor(wilson_phase(ib));
    }

    return wilson_phase;
}