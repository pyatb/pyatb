#ifndef BASE_DATA_H
#define BASE_DATA_H

#include <iostream>
#include <array>
#include "tools.h"
#include "cell_atom.h"
#include "constants.h"
#include "omp.h"

class base_data
{
    // May need to be modified
    friend class xr_operation;
    friend class interface_python;

public:
    base_data();
    ~base_data();

    int get_basis_num(){return this->basis_num;}
    double get_lattice_constant(){return this->lattice_constant;}
    double get_primitiveCell_volume(){return this->lattice_vector.determinant() * this->lattice_constant * this-> lattice_constant * this->lattice_constant;}
    Matrix3d get_lattice_vector(){return this->lattice_vector;}
    Matrix3d get_reciprocal_vector(){return this->reciprocal_vector;}

    VectorXcd get_exp_ikR_1k(const VectorXd &k_direct_coor);

    MatrixXcd get_exp_ikR(const MatrixXd &k_direct_coor);

    MatrixXcd get_Xk_triu(
        const VectorXcd &exp_ikR, 
        const MatrixXcd &XR_upperTriangleOfDenseMatrix
    );

    MatrixXcd get_Xk_triu_sparse(
        const VectorXcd &exp_ikR, 
        const SparseMatrixXcdC &XR_upperTriangleOfSparseMatrix
    );

    MatrixXcd get_partial_Xk_triu(
        const VectorXcd &exp_ikR, 
        const MatrixXcd &XR_upperTriangleOfDenseMatrix, 
        const int &partial_alpha
    );

    MatrixXcd get_partial_Xk_triu_sparse(
        const VectorXcd &exp_ikR, 
        const SparseMatrixXcdC &XR_upperTriangleOfSparseMatrix, 
        const int &partial_alpha
    );


    MatrixXcd get_partial2_Xk_triu(
        const VectorXcd &exp_ikR,
        const MatrixXcd &XR_upperTriangleOfDenseMatrix,
        const int &partial_alpha,
        const int &partial_beta
    );

    MatrixXcd get_partial2_Xk_triu_sparse(
        const VectorXcd &exp_ikR,
        const SparseMatrixXcdC &XR_upperTriangleOfSparseMatrix,
        const int &partial_alpha,
        const int &partial_beta
    );

    double lattice_constant;
    Matrix3d lattice_vector;
    Matrix3d reciprocal_vector;
    std::vector<cell_atom> atom;

private:
    int basis_num;
    int R_num;
    MatrixXd R_direct_coor;
    MatrixXd R_cartesian_coor;

    // HR, SR, rR is sparse matrix
    bool use_XR_sparse = false;

    // Dense Matrix
    MatrixXcd HR_upperTriangleOfDenseMatrix;
    MatrixXcd SR_upperTriangleOfDenseMatrix;
    std::array<MatrixXcd, 3> rR_upperTriangleOfDenseMatrix;
    
    // Sparse Matrix
    SparseMatrixXcdC HR_upperTriangleOfSparseMatrix;
    SparseMatrixXcdC SR_upperTriangleOfSparseMatrix;
    std::array<SparseMatrixXcdC, 3> rR_upperTriangleOfSparseMatrix;

    MatrixXd basis_center_position_direct;
    MatrixXd basis_center_position_cartesian;
};

#endif