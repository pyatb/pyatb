#include "tools.h"

MatrixXcd tools::convert_tril(
    const int &basis_num, 
    const VectorXcd &upperTriangleOfDenseMatrix
)
{
    MatrixXcd container;
    container.setZero(basis_num, basis_num);

    int count = 0;
    for (int row = 0; row < basis_num; row++)
    {
        for (int col = row; col < basis_num; col++)
        {
            container(row, col) = upperTriangleOfDenseMatrix(count);
            count++;
        }
    }

    return container;
}

void tools::diagonalize_SelfAdjointMatrix_eigenvaluesOnly(
    const MatrixXcd &H_k,
    VectorXd &eigenvalues
)
{
    SelfAdjointEigenSolver<MatrixXcd> eigenSolver;
    eigenSolver.compute(H_k, EigenvaluesOnly);
    eigenvalues = eigenSolver.eigenvalues();
}

void tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_1k(
    const MatrixXcd &H_k, 
    const MatrixXcd &S_k, 
    VectorXd &eigenvalues
)
{
    GeneralizedSelfAdjointEigenSolver<MatrixXcd> eigenSolver;
    eigenSolver.compute(H_k, S_k, EigenvaluesOnly);
    eigenvalues = eigenSolver.eigenvalues();
}

void tools::diagonalize_GeneralizedSelfAdjointMatrix_1k(
    const MatrixXcd &H_k, 
    const MatrixXcd &S_k, 
    MatrixXcd &eigenvectors,
    VectorXd &eigenvalues
)
{
    GeneralizedSelfAdjointEigenSolver<MatrixXcd> eigenSolver;  
    eigenSolver.compute(H_k, S_k);
    eigenvectors = eigenSolver.eigenvectors();
    eigenvalues = eigenSolver.eigenvalues();
}


void tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly(
    const std::vector<MatrixXcd> &H_k, 
    const std::vector<MatrixXcd> &S_k, 
    std::vector<VectorXd> &eigenvalues
)
{
    int kpoint_num = H_k.size();
    int max_num_threads = omp_get_max_threads();

    eigenvalues.resize(kpoint_num);

    if (max_num_threads > kpoint_num)
    {
        GeneralizedSelfAdjointEigenSolver<MatrixXcd> eigenSolver;
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            eigenSolver.compute(H_k[ik], S_k[ik], EigenvaluesOnly);
            eigenvalues[ik] = eigenSolver.eigenvalues();
        }
    }
    else
    {
        GeneralizedSelfAdjointEigenSolver<MatrixXcd> eigenSolver;

        #pragma omp parallel for private(eigenSolver) schedule(static)
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            
            eigenSolver.compute(H_k[ik], S_k[ik], EigenvaluesOnly);
            eigenvalues[ik] = eigenSolver.eigenvalues();
        }
    }
}

void tools::diagonalize_GeneralizedSelfAdjointMatrix(
    const std::vector<MatrixXcd> &H_k, 
    const std::vector<MatrixXcd> &S_k, 
    std::vector<MatrixXcd> &eigenvectors,
    std::vector<VectorXd> &eigenvalues
)
{
    int kpoint_num = H_k.size();
    int max_num_threads = omp_get_max_threads();

    eigenvectors.resize(kpoint_num);
    eigenvalues.resize(kpoint_num);

    if (max_num_threads > kpoint_num)
    {
        GeneralizedSelfAdjointEigenSolver<MatrixXcd> eigenSolver;
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            eigenSolver.compute(H_k[ik], S_k[ik]);
            eigenvectors[ik] = eigenSolver.eigenvectors();
            eigenvalues[ik] = eigenSolver.eigenvalues();
        }
    }
    else
    {
        GeneralizedSelfAdjointEigenSolver<MatrixXcd> eigenSolver;

        #pragma omp parallel for private(eigenSolver) schedule(static)
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            
            eigenSolver.compute(H_k[ik], S_k[ik]);
            eigenvectors[ik] = eigenSolver.eigenvectors();
            eigenvalues[ik] = eigenSolver.eigenvalues();
        }
    }
}