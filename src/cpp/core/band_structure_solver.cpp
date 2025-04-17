#include "band_structure_solver.h"


void band_structure_solver::get_eigenvalues_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    VectorXd &eigenvalues
)
{
    MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR);
    MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR);
    tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_1k(Hk, Sk, eigenvalues);
}

void band_structure_solver::get_eigenvalues_range_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &lower_band_index, // counting from 1.
    const int &upper_band_index, // counting from 1, upper_band_index >= lower_band_index
    VectorXd &eigenvalues
)
{
    MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR);
    MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR);
    tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_range_1k(Hk, Sk, lower_band_index, upper_band_index, eigenvalues);
}


void band_structure_solver::get_eigenvalues_eigenvectors_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    VectorXd &eigenvalues, 
    MatrixXcd &eigenvectors
)
{
    MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR);
    MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR);
    tools::diagonalize_GeneralizedSelfAdjointMatrix_1k(Hk, Sk, eigenvectors, eigenvalues);
}


void band_structure_solver::get_eigenvalues_eigenvectors_range_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &lower_band_index, // counting from 1.
    const int &upper_band_index, // counting from 1, upper_band_index >= lower_band_index
    VectorXd &eigenvalues, 
    MatrixXcd &eigenvectors
)
{
    MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR);
    MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR);
    tools::diagonalize_GeneralizedSelfAdjointMatrix_range_1k(Hk, Sk, lower_band_index, upper_band_index, eigenvectors, eigenvalues);
}


void band_structure_solver::get_eigenvalues(
    base_data &Base_Data, 
    const MatrixXd &k_direct_coor, 
    std::vector<VectorXd> &eigenvalues
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR.row(ik));
        MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR.row(ik));
        tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_1k(Hk, Sk, eigenvalues[ik]);
    }
}


void band_structure_solver::get_eigenvalues_range(
    base_data &Base_Data, 
    const MatrixXd &k_direct_coor, 
    const int &lower_band_index, // counting from 1.
    const int &upper_band_index, // counting from 1, upper_band_index >= lower_band_index
    std::vector<VectorXd> &eigenvalues
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR.row(ik));
        MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR.row(ik));
        tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_range_1k(Hk, Sk, lower_band_index, upper_band_index, eigenvalues[ik]);
    }
}


void band_structure_solver::get_eigenvalues_eigenvectors(
    base_data &Base_Data, 
    const MatrixXd &k_direct_coor, 
    std::vector<VectorXd> &eigenvalues, 
    std::vector<MatrixXcd> &eigenvectors
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR.row(ik));
        MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR.row(ik));
        tools::diagonalize_GeneralizedSelfAdjointMatrix_1k(Hk, Sk, eigenvectors[ik], eigenvalues[ik]);
    }
}


void band_structure_solver::get_eigenvalues_eigenvectors_range(
    base_data &Base_Data, 
    const MatrixXd &k_direct_coor, 
    const int &lower_band_index, // counting from 1.
    const int &upper_band_index, // counting from 1, upper_band_index >= lower_band_index
    std::vector<VectorXd> &eigenvalues, 
    std::vector<MatrixXcd> &eigenvectors
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR.row(ik));
        MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR.row(ik));
        tools::diagonalize_GeneralizedSelfAdjointMatrix_range_1k(Hk, Sk, lower_band_index, upper_band_index, eigenvectors[ik], eigenvalues[ik]);
    }
}