#include "velocity_solver.h"

MatrixXcd velocity_solver::cal_velocity_1k_base(
    base_data &Base_Data,
    const VectorXcd &exp_ikR, 
    const VectorXd &eigenvalues, 
    const MatrixXcd &eigenvectors, 
    const int &alpha
)
{
    MatrixXcd partial_Hk = eigenvectors.adjoint() * xr_operation::get_partial_Hk(Base_Data, exp_ikR, alpha) * eigenvectors;
    MatrixXcd partial_Sk = eigenvectors.adjoint() * xr_operation::get_partial_Sk(Base_Data, exp_ikR, alpha) * eigenvectors;
    MatrixXcd rk = eigenvectors.adjoint() * xr_operation::get_rk(Base_Data, exp_ikR, alpha) * eigenvectors;

    int basis_num = Base_Data.get_basis_num();
    MatrixXcd velocity_matrix = MatrixXcd::Zero(basis_num, basis_num);

    for (int row = 0; row < basis_num; row++)
    {
        for (int col = row; col < basis_num; col++)
        {
            velocity_matrix(row, col) =  partial_Hk(row, col) - eigenvalues(row) * partial_Sk(row, col)
                                       + IMAG_UNIT * (eigenvalues(row) - eigenvalues(col)) * rk(row, col);
        }
    }
    
    velocity_matrix = velocity_matrix.selfadjointView<Upper>();
    
    return velocity_matrix;
}


void velocity_solver::get_velocity_matrix_alpha(
    base_data &Base_Data,
    const MatrixXd &k_direct_coor, 
    const int &alpha, 
    std::vector<MatrixXcd> &velocity_matrix
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);
    velocity_matrix.resize(kpoint_num);

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXd eigenvalues;
        MatrixXcd eigenvectors;
        band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors);
        velocity_matrix[ik] = velocity_solver::cal_velocity_1k_base(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors, alpha);
    }

}

void velocity_solver::get_velocity_matrix(
    base_data &Base_Data,
    const MatrixXd &k_direct_coor, 
    std::array<std::vector<MatrixXcd>, 3> &velocity_matrix
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);
    for (int i = 0; i < 3; ++i)
    {
        velocity_matrix[i].resize(kpoint_num);
    }

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXd eigenvalues;
        MatrixXcd eigenvectors;
        band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors);

        for (int i = 0; i < 3; ++i)
        {
            velocity_matrix[i][ik] = velocity_solver::cal_velocity_1k_base(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors, i);
        }
        
    }
}
