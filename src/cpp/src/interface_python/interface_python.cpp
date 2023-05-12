#include "interface_python.h"
#include "../core/xr_operation.h"
#include "../core/tools.h"
#include "../core/berry_curvature_solver.h"
#include "../core/berry_phase_solver.h"
#include "../core/optical_conductivity_solver.h"
#include "../core/velocity_solver.h"
#include "../core/cell_atom.h"
#include "../core/bandunfolding_solver.h"
#include "../core/shift_current_solver.h"

interface_python::interface_python(double lattice_constant, Matrix3d &lattice_vector)
{
    Base_Data.lattice_constant = lattice_constant;
    Base_Data.lattice_vector.swap(lattice_vector);
    Base_Data.reciprocal_vector = (Base_Data.lattice_vector.inverse()).transpose();
}

interface_python::~interface_python(){}

void interface_python::set_HSR(
    int R_num,
    MatrixXd &R_direct_coor,
    int basis_num,
    MatrixXcd &HR_upperTriangleOfDenseMatrix,
    MatrixXcd &SR_upperTriangleOfDenseMatrix
)
{
    // R_num, R_direct_coor
    Base_Data.R_num = R_num;
    Base_Data.R_direct_coor.swap(R_direct_coor);

    // calculate R_cartesian_coor from R_direct_coor
    Base_Data.R_cartesian_coor = Base_Data.R_direct_coor * Base_Data.lattice_vector * Base_Data.lattice_constant;
    
    // basis num
    Base_Data.basis_num = basis_num;

    // HR, SR
    Base_Data.HR_upperTriangleOfDenseMatrix.swap(HR_upperTriangleOfDenseMatrix);
    Base_Data.SR_upperTriangleOfDenseMatrix.swap(SR_upperTriangleOfDenseMatrix);
}


void interface_python::set_HSR_sparse(
    int R_num, 
    MatrixXd &R_direct_coor, 
    int basis_num,
    SparseMatrixXcdC &HR_upperTriangleOfSparseMatrix, 
    SparseMatrixXcdC &SR_upperTriangleOfSparseMatrix
)
{
    // R_num, R_direct_coor
    Base_Data.R_num = R_num;
    Base_Data.R_direct_coor.swap(R_direct_coor);

    // calculate R_cartesian_coor from R_direct_coor
    Base_Data.R_cartesian_coor = Base_Data.R_direct_coor * Base_Data.lattice_vector * Base_Data.lattice_constant;
    
    // basis num
    Base_Data.basis_num = basis_num;

    // HR, SR
    Base_Data.HR_upperTriangleOfSparseMatrix.swap(HR_upperTriangleOfSparseMatrix);
    Base_Data.SR_upperTriangleOfSparseMatrix.swap(SR_upperTriangleOfSparseMatrix);

    // set sparse matrix flag
    Base_Data.use_XR_sparse = true;
}


void interface_python::set_rR(
    MatrixXcd &rR_x,
    MatrixXcd &rR_y,
    MatrixXcd &rR_z
)
{
    Base_Data.rR_upperTriangleOfDenseMatrix[0].swap(rR_x);
    Base_Data.rR_upperTriangleOfDenseMatrix[1].swap(rR_y);
    Base_Data.rR_upperTriangleOfDenseMatrix[2].swap(rR_z);

    Base_Data.basis_center_position_cartesian.setZero(Base_Data.basis_num, 3);

    Vector3d R_zero(0.0, 0.0, 0.0);
    vector<int> diagonal_index(Base_Data.basis_num);
    diagonal_index[0] = 0;
    for (int i = 1; i < Base_Data.basis_num; ++i)
    {
        diagonal_index[i] = diagonal_index[i-1] + Base_Data.basis_num - i + 1;
    }
    for (int iR = 0; iR < Base_Data.R_num; ++iR)
    {
        if (Base_Data.R_direct_coor.row(iR).isApprox(R_zero.transpose()))
        {
            for (int row = 0; row < Base_Data.basis_num; ++row)
            {
                Base_Data.basis_center_position_cartesian(row, 0) = Base_Data.rR_upperTriangleOfDenseMatrix[0](iR, diagonal_index[row]).real();
                Base_Data.basis_center_position_cartesian(row, 1) = Base_Data.rR_upperTriangleOfDenseMatrix[1](iR, diagonal_index[row]).real();
                Base_Data.basis_center_position_cartesian(row, 2) = Base_Data.rR_upperTriangleOfDenseMatrix[2](iR, diagonal_index[row]).real();
            }

            break;
        }
    }

    Base_Data.basis_center_position_direct = Base_Data.basis_center_position_cartesian * Base_Data.lattice_vector.inverse() / Base_Data.lattice_constant;

}

void interface_python::set_rR_sparse(
    SparseMatrixXcdC &rR_x,
    SparseMatrixXcdC &rR_y,
    SparseMatrixXcdC &rR_z
)
{
    Base_Data.rR_upperTriangleOfSparseMatrix[0].swap(rR_x);
    Base_Data.rR_upperTriangleOfSparseMatrix[1].swap(rR_y);
    Base_Data.rR_upperTriangleOfSparseMatrix[2].swap(rR_z);

    Base_Data.basis_center_position_cartesian.setZero(Base_Data.basis_num, 3);

    Vector3d R_zero(0.0, 0.0, 0.0);
    vector<int> diagonal_index(Base_Data.basis_num);
    diagonal_index[0] = 0;
    for (int i = 1; i < Base_Data.basis_num; ++i)
    {
        diagonal_index[i] = diagonal_index[i-1] + Base_Data.basis_num - i + 1;
    }
    for (int iR = 0; iR < Base_Data.R_num; ++iR)
    {
        if (Base_Data.R_direct_coor.row(iR).isApprox(R_zero.transpose()))
        {
            for (int row = 0; row < Base_Data.basis_num; ++row)
            {
                Base_Data.basis_center_position_cartesian(row, 0) = Base_Data.rR_upperTriangleOfSparseMatrix[0].coeff(iR, diagonal_index[row]).real();
                Base_Data.basis_center_position_cartesian(row, 1) = Base_Data.rR_upperTriangleOfSparseMatrix[1].coeff(iR, diagonal_index[row]).real();
                Base_Data.basis_center_position_cartesian(row, 2) = Base_Data.rR_upperTriangleOfSparseMatrix[2].coeff(iR, diagonal_index[row]).real();
            }

            break;
        }
    }

    Base_Data.basis_center_position_direct = Base_Data.basis_center_position_cartesian * Base_Data.lattice_vector.inverse() / Base_Data.lattice_constant;

}

void interface_python::set_single_atom_position(
    std::string atom_label,
    int na,
    MatrixXd &tau_car
)
{
    Base_Data.atom.push_back(cell_atom());
    cell_atom &temp_atom = Base_Data.atom.back();
    temp_atom.set_atom_position(
        Base_Data.lattice_constant,
        Base_Data.lattice_vector,
        atom_label,
        na,
        tau_car
    );
}

void interface_python::set_single_atom_orb(
    int &atom_index,
    int &nwl,
    std::vector<int> &l_nchi,
    int &mesh,
    double &dr,
    MatrixXd &numerical_orb
)
{
    cell_atom &temp_atom = Base_Data.atom[atom_index];
    temp_atom.set_numerical_orb(
        nwl,
        l_nchi,
        mesh,
        dr,
        numerical_orb
    );
}


MatrixXcd& interface_python::get_HR()
{
    return Base_Data.HR_upperTriangleOfDenseMatrix;
}


MatrixXcd& interface_python::get_SR()
{
    return Base_Data.SR_upperTriangleOfDenseMatrix;
}


SparseMatrixXcdC& interface_python::get_HR_sparse()
{
    return Base_Data.HR_upperTriangleOfSparseMatrix;
}


SparseMatrixXcdC& interface_python::get_SR_sparse()
{
    return Base_Data.SR_upperTriangleOfSparseMatrix;
}


MatrixXcd& interface_python::get_rR(int direction)
{
    return Base_Data.rR_upperTriangleOfDenseMatrix[direction];
}


SparseMatrixXcdC& interface_python::get_rR_sparse(int direction)
{
    return Base_Data.rR_upperTriangleOfSparseMatrix[direction];
}


void interface_python::update_HR_sparse(SparseMatrixXcdC &HR)
{
    Base_Data.HR_upperTriangleOfSparseMatrix.swap(HR);
}


void interface_python::update_SR_sparse(SparseMatrixXcdC &SR)
{
    Base_Data.SR_upperTriangleOfSparseMatrix.swap(SR);
}


void interface_python::update_rR_sparse(int direction, SparseMatrixXcdC &rR_d)
{
    Base_Data.rR_upperTriangleOfSparseMatrix[direction].swap(rR_d);
}


void interface_python::get_Hk(
    const MatrixXd &k_direct_coor, 
    py::array_t<std::complex<double>> &Hk
)
{
    auto Hk_data = Hk.mutable_unchecked<3>();

    const int kpoint_num = k_direct_coor.rows();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);
    int max_num_threads = omp_get_max_threads();

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd temp_Hk = xr_operation::get_Hk(Base_Data, exp_ikR.row(ik));

        for (int row = 0; row < Base_Data.basis_num; ++row)
        {
            for (int col = 0; col < Base_Data.basis_num; ++col)
            {
                Hk_data(ik, row, col) = temp_Hk(row, col);
            }
        }
    }

}


void interface_python::get_Sk(
    const MatrixXd &k_direct_coor, 
    py::array_t<std::complex<double>> &Sk
)
{
    auto Sk_data = Sk.mutable_unchecked<3>();

    const int kpoint_num = k_direct_coor.rows();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);
    int max_num_threads = omp_get_max_threads();

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        MatrixXcd temp_Sk = xr_operation::get_Sk(Base_Data, exp_ikR.row(ik));

        for (int row = 0; row < Base_Data.basis_num; ++row)
        {
            for (int col = 0; col < Base_Data.basis_num; ++col)
            {
                Sk_data(ik, row, col) = temp_Sk(row, col);
            }
        }
    }

}


// void interface_python::get_HSk_surface(
//     int direction, 
//     int coupling_layers,
//     const MatrixXd &k_direct_coor, 
//     py::array_t<std::complex<double>> &Hk00, 
//     py::array_t<std::complex<double>> &Hk01, 
//     py::array_t<std::complex<double>> &Sk00, 
//     py::array_t<std::complex<double>> &Sk01
// )
// {
//     auto Hk00_data = Hk00.mutable_unchecked<3>();
//     auto Hk01_data = Hk01.mutable_unchecked<3>();
//     auto Sk00_data = Sk00.mutable_unchecked<3>();
//     auto Sk01_data = Sk01.mutable_unchecked<3>();

//     const int kpoint_num = k_direct_coor.rows();
//     vector<MatrixXcd> tem_Hk00(kpoint_num);
//     vector<MatrixXcd> tem_Hk01(kpoint_num);
//     vector<MatrixXcd> tem_Sk00(kpoint_num);
//     vector<MatrixXcd> tem_Sk01(kpoint_num);

//     int matrix_dim = coupling_layers * Base_Data.basis_num;
//     for (int ik = 0; ik < kpoint_num; ++ik)
//     {
//         tem_Hk00[ik].setZero(matrix_dim, matrix_dim);
//         tem_Hk01[ik].setZero(matrix_dim, matrix_dim);
//         tem_Sk00[ik].setZero(matrix_dim, matrix_dim);
//         tem_Sk01[ik].setZero(matrix_dim, matrix_dim);
//     }

//     Base_Data.SurfaceState_Xk(direction, k_direct_coor, coupling_layers, 'H', tem_Hk00, tem_Hk01);
//     Base_Data.SurfaceState_Xk(direction, k_direct_coor, coupling_layers, 'S', tem_Sk00, tem_Sk01);

//     for (int ik = 0; ik < kpoint_num; ++ik)
//     {
//         for (int row = 0; row < matrix_dim; ++row)
//         {
//             for (int col = 0; col < matrix_dim; ++col)
//             {
//                 Hk00_data(ik, row, col) = tem_Hk00[ik](row, col);
//                 Hk01_data(ik, row, col) = tem_Hk01[ik](row, col);
//                 Sk00_data(ik, row, col) = tem_Sk00[ik](row, col);
//                 Sk01_data(ik, row, col) = tem_Sk01[ik](row, col);
//             }
//         }
//     }

// }


void interface_python::diago_H(
    const MatrixXd &k_direct_coor,
    py::array_t<std::complex<double>> &eigenvectors,
    py::array_t<double> &eigenvalues
)
{
    auto eigenvectors_data = eigenvectors.mutable_unchecked<3>();
    auto eigenvalues_data = eigenvalues.mutable_unchecked<2>();

    const int kpoint_num = k_direct_coor.rows();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);
    int max_num_threads = omp_get_max_threads();

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXd temp_eigenvalues;
        MatrixXcd temp_eigenvectors;
        band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), temp_eigenvalues, temp_eigenvectors);

        for (int row = 0; row < Base_Data.basis_num; ++row)
        {
            for (int col = 0; col < Base_Data.basis_num; ++col)
            {
                eigenvectors_data(ik, row, col) = temp_eigenvectors(row, col);
            }

            eigenvalues_data(ik, row) = temp_eigenvalues[row];
        }
    }

}

void interface_python::diago_H_eigenvaluesOnly(
    const MatrixXd &k_direct_coor,
    py::array_t<double> &eigenvalues
)
{
    auto eigenvalues_data = eigenvalues.mutable_unchecked<2>();
    const int kpoint_num = k_direct_coor.rows();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);
    int max_num_threads = omp_get_max_threads();

    #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXd temp_eigenvalues;
        band_structure_solver::get_eigenvalues_1k(Base_Data, exp_ikR.row(ik), temp_eigenvalues);

        for (int row = 0; row < Base_Data.basis_num; ++row)
        {
            eigenvalues_data(ik, row) = temp_eigenvalues[row];
        }
    }
}


void interface_python::get_total_berry_curvature_fermi(
    const MatrixXd &k_direct_coor,
    const double &fermi_energy,
    const int mode,
    py::array_t<double> &total_berry_curvature
)
{
    auto total_berry_curvature_data = total_berry_curvature.mutable_unchecked<2>();
    const int kpoint_num = k_direct_coor.rows();
    MatrixXd tem_berry_curvature_values = berry_curvature_solver::get_total_bc_fermi(Base_Data, k_direct_coor, fermi_energy, mode);

    for (int ik = 0; ik < kpoint_num; ++ik)
    {       
        total_berry_curvature_data(ik, 0) = tem_berry_curvature_values(ik, 0);
        total_berry_curvature_data(ik, 1) = tem_berry_curvature_values(ik, 1);
        total_berry_curvature_data(ik, 2) = tem_berry_curvature_values(ik, 2);
    }
}

void interface_python::get_total_berry_curvature_occupiedNumber(
    const MatrixXd &k_direct_coor,
    const int &occupied_band_num,
    const int mode,
    py::array_t<double> &total_berry_curvature
)
{
    auto total_berry_curvature_data = total_berry_curvature.mutable_unchecked<2>();
    const int kpoint_num = k_direct_coor.rows();
    MatrixXd tem_berry_curvature_values = berry_curvature_solver::get_total_bc_occupiedNumber(Base_Data, k_direct_coor, occupied_band_num, mode);

    for (int ik = 0; ik < kpoint_num; ++ik)
    {       
        total_berry_curvature_data(ik, 0) = tem_berry_curvature_values(ik, 0);
        total_berry_curvature_data(ik, 1) = tem_berry_curvature_values(ik, 1);
        total_berry_curvature_data(ik, 2) = tem_berry_curvature_values(ik, 2);
    }
}


// void interface_python::get_berry_curvature_and_eigenvalues_by_fermi(
//     const MatrixXd &k_direct_coor,
//     py::array_t<double> &berry_curvature_values, 
//     py::array_t<double> &eigenvalues,
//     const double &fermi_energy,
//     const int mode
// )
// {
//     auto berry_curvature_values_data = berry_curvature_values.mutable_unchecked<2>();
//     auto eigenvalues_data = eigenvalues.mutable_unchecked<2>();

//     const int kpoint_num = k_direct_coor.rows();
//     vector<array<double, 3>> tem_berry_curvature_values(kpoint_num);
//     berry_curvature_solver BCSolver(Base_Data, k_direct_coor, fermi_energy);
//     BCSolver.get_berry_curvature(tem_berry_curvature_values, mode);
//     const vector<VectorXd> &tem_eigenvalues = BCSolver.get_eigenvalues();

//     for (int ik = 0; ik < kpoint_num; ++ik)
//     {       
//         berry_curvature_values_data(ik, 0) = tem_berry_curvature_values[ik][0];
//         berry_curvature_values_data(ik, 1) = tem_berry_curvature_values[ik][1];
//         berry_curvature_values_data(ik, 2) = tem_berry_curvature_values[ik][2];

//         for (int row = 0; row < Base_Data.basis_num; ++row)
//         {
//             eigenvalues_data(ik, row) = tem_eigenvalues[ik][row];
//         }
//     }

// }

// void interface_python::get_berry_curvature_and_eigenvalues_by_occupy(
//     const MatrixXd &k_direct_coor,
//     py::array_t<double> &berry_curvature_values, 
//     py::array_t<double> &eigenvalues,
//     const int &occupied_band_num,
//     const int mode
// )
// {
//     auto berry_curvature_values_data = berry_curvature_values.mutable_unchecked<2>();
//     auto eigenvalues_data = eigenvalues.mutable_unchecked<2>();

//     const int kpoint_num = k_direct_coor.rows();
//     vector<array<double, 3>> tem_berry_curvature_values(kpoint_num);
//     berry_curvature_solver BCSolver(Base_Data, k_direct_coor, occupied_band_num);
//     BCSolver.get_berry_curvature(tem_berry_curvature_values, mode);
//     const vector<VectorXd> &tem_eigenvalues = BCSolver.get_eigenvalues();

//     for (int ik = 0; ik < kpoint_num; ++ik)
//     {       
//         berry_curvature_values_data(ik, 0) = tem_berry_curvature_values[ik][0];
//         berry_curvature_values_data(ik, 1) = tem_berry_curvature_values[ik][1];
//         berry_curvature_values_data(ik, 2) = tem_berry_curvature_values[ik][2];

//         for (int row = 0; row < Base_Data.basis_num; ++row)
//         {
//             eigenvalues_data(ik, row) = tem_eigenvalues[ik][row];
//         }
//     }

// }


double interface_python::get_berry_phase_of_loop(
    const MatrixXd &k_direct_coor_loop, 
    const int &occupied_band_num
)
{
    double phase = berry_phase_solver::get_berry_phase_of_loop(Base_Data, k_direct_coor_loop, occupied_band_num);
    return phase;
}


VectorXd interface_python::get_wilson_loop(
    const MatrixXd &k_direct_coor_loop, 
    const int &occupied_band_num
)
{
    VectorXd tem_wilson_phase = berry_phase_solver::calculate_wilson_loop(Base_Data, k_direct_coor_loop, occupied_band_num);
    return tem_wilson_phase;
}


void interface_python::get_optical_conductivity_by_kubo(
    const int &nspin,
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const double &eta,
    const int &occupied_band_num,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num,
    const int &method,
    py::array_t<std::complex<double>> optical_conductivity,
    py::array_t<std::complex<double>> dielectric_function
)
{
    auto oc_data = optical_conductivity.mutable_unchecked<2>();
    auto df_data = dielectric_function.mutable_unchecked<2>();

    MatrixXcd oc_tem(9, omega_num);
    oc_tem.setZero();

    MatrixXcd df_tem(9, omega_num);
    df_tem.setZero();

    optical_conductivity_solver OCSolver;
    OCSolver.set_parameters(nspin, omega_num, domega, start_omega, eta, occupied_band_num, k_direct_coor, total_kpoint_num);
    OCSolver.get_optical_conductivity_by_kubo(Base_Data, method, oc_tem, df_tem);

    for (int i = 0; i < 9; ++i)
    {
        for (int i_omega = 0; i_omega < omega_num; ++i_omega)
        {
            oc_data(i, i_omega) += oc_tem(i, i_omega);
            df_data(i, i_omega) += df_tem(i, i_omega);
        }
    }

}

void interface_python::get_shift_current(
    const int &nspin,
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const int &smearing_method,
    const double &eta,
    const int &occupied_band_num,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num,
    const int &method,
    py::array_t<double> shift_current
)
{
    auto data = shift_current.mutable_unchecked<2>();

    shift_current_solver SCS;
    SCS.set_parameters(nspin, omega_num, domega, start_omega, smearing_method, eta);
    MatrixXd tem = SCS.get_shift_current_conductivity(Base_Data, k_direct_coor, total_kpoint_num, occupied_band_num, method);

    for (int i = 0; i < 18; ++i)
    {
        for (int i_omega = 0; i_omega < omega_num; ++i_omega)
        {
            data(i, i_omega) += tem(i, i_omega);
        }
    }
}

void interface_python::get_velocity_matrix(
    const MatrixXd &k_direct_coor,
    py::array_t<double> &eigenvalues,
    py::array_t<std::complex<double>> &velocity_matrix
)
{
    auto eigenvalues_data = eigenvalues.mutable_unchecked<2>();
    auto velocity_matrix_data = velocity_matrix.mutable_unchecked<4>();
    const int kpoint_num = k_direct_coor.rows();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);

    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXd temp_eigenvalues;
        MatrixXcd temp_eigenvectors;
        band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), temp_eigenvalues, temp_eigenvectors);

        MatrixXcd temp_velocity[3];
        for (int a = 0; a < 3; ++a)
        {
            temp_velocity[a] = velocity_solver::cal_velocity_1k_base(Base_Data, exp_ikR.row(ik), temp_eigenvalues, temp_eigenvectors, a);
        }

        for (int row = 0; row < Base_Data.basis_num; ++row)
        {
            eigenvalues_data(ik, row) = temp_eigenvalues(row);
        }

        for (int a = 0; a < 3; ++a)
        {
            for (int row = 0; row < Base_Data.basis_num; ++row)
            {
                for (int col = 0; col < Base_Data.basis_num; ++col)
                {
                    velocity_matrix_data(ik, a, row, col) = temp_velocity[a](row, col);
                }
            }
        }
    }

}


void interface_python::get_bandunfolding(
    const Matrix3d &M_matrix,
    const MatrixXd &kvect_direct,
    const double &ecut,
    const int &min_bandindex,
    const int &max_bandindex,
    const int &nspin,
    py::array_t<double> &P,
    py::array_t<double> &E
)
{
    auto P_data = P.mutable_unchecked<2>();
    auto E_data = E.mutable_unchecked<2>();
    const int kpoint_num = kvect_direct.rows();
    const int select_band_num = max_bandindex - min_bandindex + 1;

    bandunfolding_solver BF_solver;
    BF_solver.set_M_matrix(Base_Data.lattice_constant, Base_Data.lattice_vector, M_matrix);
    MatrixXd temp_P, temp_E;
    BF_solver.output_spectral_weight(Base_Data, kvect_direct, ecut, min_bandindex, max_bandindex, nspin, temp_P, temp_E);

    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        for (int ib = 0; ib < select_band_num; ++ib)
        {
            P_data(ik, ib) = temp_P(ik, ib);
            E_data(ik, ib) = temp_E(ik, ib);
        }
    }
}


PYBIND11_MODULE(interface_python, m)
{
    py::class_<interface_python>(m, "interface_python")
        .def(py::init<double, Matrix3d &>())
        // .def("init", &interface_python::init)
        .def("set_HSR", &interface_python::set_HSR)
        .def("set_HSR_sparse", &interface_python::set_HSR_sparse)
        .def("set_rR", &interface_python::set_rR)
        .def("set_rR_sparse", &interface_python::set_rR_sparse)
        .def("set_single_atom_position", &interface_python::set_single_atom_position)
        .def("set_single_atom_orb", &interface_python::set_single_atom_orb)
        .def("get_HR", &interface_python::get_HR, py::return_value_policy::reference_internal)
        .def("get_SR", &interface_python::get_SR, py::return_value_policy::reference_internal)
        .def("get_HR_sparse", &interface_python::get_HR_sparse, py::return_value_policy::reference_internal)
        .def("get_SR_sparse", &interface_python::get_SR_sparse, py::return_value_policy::reference_internal)
        .def("get_rR", &interface_python::get_rR, py::return_value_policy::reference_internal)
        .def("get_rR_sparse", &interface_python::get_rR_sparse, py::return_value_policy::reference_internal)
        .def("update_HR_sparse", &interface_python::update_HR_sparse)
        .def("update_SR_sparse", &interface_python::update_SR_sparse)
        .def("update_rR_sparse", &interface_python::update_rR_sparse)
        .def("get_Hk", &interface_python::get_Hk)
        .def("get_Sk", &interface_python::get_Sk)
        // .def("get_HSk_surface", &interface_python::get_HSk_surface)
        .def("diago_H", &interface_python::diago_H)
        .def("diago_H_eigenvaluesOnly", &interface_python::diago_H_eigenvaluesOnly)
        .def("get_total_berry_curvature_fermi", &interface_python::get_total_berry_curvature_fermi)
        .def("get_total_berry_curvature_occupiedNumber", &interface_python::get_total_berry_curvature_occupiedNumber)
        .def("get_berry_phase_of_loop", &interface_python::get_berry_phase_of_loop)
        .def("get_wilson_loop", &interface_python::get_wilson_loop)
        .def("get_optical_conductivity_by_kubo", &interface_python::get_optical_conductivity_by_kubo)
        .def("get_shift_current", &interface_python::get_shift_current)
        .def("get_velocity_matrix", &interface_python::get_velocity_matrix)
        .def("get_bandunfolding", &interface_python::get_bandunfolding);
}