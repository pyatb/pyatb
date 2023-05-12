#ifndef TOOL_H
#define TOOL_H

#include <vector>
#include "use_eigen.h"
#include "omp.h"

using namespace Eigen;

class tools
{
public:
    static MatrixXcd convert_tril(
        const int &basis_num, 
        const VectorXcd &upperTriangleOfDenseMatrix
    );

    static void diagonalize_SelfAdjointMatrix_eigenvaluesOnly(
        const MatrixXcd &H_k,
        VectorXd &eigenvalues
    );

    static void diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_1k(
        const MatrixXcd &H_k, 
        const MatrixXcd &S_k, 
        VectorXd &eigenvalues
    );

    static void diagonalize_GeneralizedSelfAdjointMatrix_1k(
        const MatrixXcd &H_k, 
        const MatrixXcd &S_k, 
        MatrixXcd &eigenvectors,
        VectorXd &eigenvalues
    );

    static void diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly(
        const std::vector<MatrixXcd> &H_k, 
        const std::vector<MatrixXcd> &S_k, 
        std::vector<VectorXd> &eigenvalues
    );

    static void diagonalize_GeneralizedSelfAdjointMatrix(
        const std::vector<MatrixXcd> &H_k, 
        const std::vector<MatrixXcd> &S_k, 
        std::vector<MatrixXcd> &eigenvectors,
        std::vector<VectorXd> &eigenvalues
    );
};



#endif