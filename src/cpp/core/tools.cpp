#include "tools.h"
#include <iostream>

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

void tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_range_1k(
    const MatrixXcd &H_k, 
    const MatrixXcd &S_k, 
    const int &lower_eigen_index, // counting from 1.
    const int &upper_eigen_index, // counting from 1, upper_band_index >= lower_band_index
    VectorXd &eigenvalues
)
{
    MatrixXcd copy_H_k = H_k;
    MatrixXcd copy_S_k = S_k;

    int itype = 1; 
    char jobz = 'N'; 
    char range = 'I'; 
    double vl = 0.0;
    double vu = 0.0;
    lapack_int il = lower_eigen_index;
    lapack_int iu = upper_eigen_index;
    double abstol = 0.0;
    lapack_int n = H_k.rows();
    lapack_int lda = n;
    lapack_int ldb = n;

    std::vector<std::complex<double>> z(n * n);
    std::vector<double> w(n);
    std::vector<lapack_int> ifail(n);

    lapack_int m;
    lapack_int info;

    info = LAPACKE_zhegvx(
        LAPACK_COL_MAJOR,
        itype,
        jobz,
        range,
        'U',
        n,
        reinterpret_cast<lapack_complex_double*>(copy_H_k.data()), lda,
        reinterpret_cast<lapack_complex_double*>(copy_S_k.data()), ldb,
        vl, vu, il, iu,
        abstol,
        &m,
        w.data(),
        reinterpret_cast<lapack_complex_double*>(z.data()), n,
        &ifail[0]
    );

    if (info > 0) {
        std::cerr << "tools.cpp: The algorithm failed to compute eigenvalues." << std::endl;
    } else if (info < 0) {
        std::cerr << "tools.cpp: Argument " << -info << " had an illegal value." << std::endl;
    }

    eigenvalues = VectorXd::Zero(m);
    for(int i = 0; i < m; ++i)
    {
        eigenvalues[i] = w[i];
    }

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

void tools::diagonalize_GeneralizedSelfAdjointMatrix_range_1k(
    const MatrixXcd &H_k, 
    const MatrixXcd &S_k, 
    const int &lower_eigen_index, // counting from 1.
    const int &upper_eigen_index, // counting from 1, upper_band_index >= lower_band_index
    MatrixXcd &eigenvectors,
    VectorXd &eigenvalues
)
{
    MatrixXcd copy_H_k = H_k;
    MatrixXcd copy_S_k = S_k;

    int itype = 1; 
    char jobz = 'V'; 
    char range = 'I'; 
    double vl = 0.0;
    double vu = 0.0;
    lapack_int il = lower_eigen_index;
    lapack_int iu = upper_eigen_index;
    double abstol = 0.0;
    lapack_int n = H_k.rows();
    lapack_int lda = n;
    lapack_int ldb = n;

    std::vector<std::complex<double>> z(n * n);
    std::vector<double> w(n);
    std::vector<lapack_int> ifail(n);

    lapack_int m;
    lapack_int info;

    info = LAPACKE_zhegvx(
        LAPACK_COL_MAJOR,
        itype,
        jobz,
        range,
        'U',
        n,
        reinterpret_cast<lapack_complex_double*>(copy_H_k.data()), lda,
        reinterpret_cast<lapack_complex_double*>(copy_S_k.data()), ldb,
        vl, vu, il, iu,
        abstol,
        &m,
        w.data(),
        reinterpret_cast<lapack_complex_double*>(z.data()), n,
        &ifail[0]
    );

    if (info > 0) {
        std::cerr << "tools.cpp: The algorithm failed to compute eigenvalues." << std::endl;
    } else if (info < 0) {
        std::cerr << "tools.cpp: Argument " << -info << " had an illegal value." << std::endl;
    }

    eigenvalues = VectorXd::Zero(m);
    for(int i = 0; i < m; ++i)
    {
        eigenvalues[i] = w[i];
    }

    eigenvectors = MatrixXcd::Zero(n, m);
    for(int i = 0; i < m; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            eigenvectors(j, i) = z[i * n + j];
        }
    }

    // Check for unconverged eigenvectors
    bool has_ifail = false;
    for(int i = 0; i < m; ++i)
    {
        if(ifail[i] > 0)
        {
            if(!has_ifail)
            {
                std::cout << "tools.cpp: The following eigenvectors failed to converge:" << std::endl;
                has_ifail = true;
            }
            std::cout << "Eigenvector " << i+1 << " failed to converge." << std::endl;
        }
    }

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

void tools::diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_range(
    const std::vector<MatrixXcd> &H_k, 
    const std::vector<MatrixXcd> &S_k, 
    const int &lower_eigen_index, // counting from 1.
    const int &upper_eigen_index, // counting from 1, upper_band_index >= lower_band_index
    std::vector<VectorXd> &eigenvalues
)
{
    int kpoint_num = H_k.size();
    int max_num_threads = omp_get_max_threads();

    eigenvalues.resize(kpoint_num);

    if (max_num_threads > kpoint_num)
    {
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_range_1k(H_k[ik], S_k[ik], lower_eigen_index, upper_eigen_index, eigenvalues[ik]);
        }
    }
    else
    {
        #pragma omp parallel for schedule(static)
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            diagonalize_GeneralizedSelfAdjointMatrix_eigenvaluesOnly_range_1k(H_k[ik], S_k[ik], lower_eigen_index, upper_eigen_index, eigenvalues[ik]);
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

void tools::diagonalize_GeneralizedSelfAdjointMatrix_range(
    const std::vector<MatrixXcd> &H_k, 
    const std::vector<MatrixXcd> &S_k, 
    const int &lower_eigen_index, // counting from 1.
    const int &upper_eigen_index, // counting from 1, upper_band_index >= lower_band_index
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
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            diagonalize_GeneralizedSelfAdjointMatrix_range_1k(H_k[ik], S_k[ik], lower_eigen_index, upper_eigen_index, eigenvectors[ik], eigenvalues[ik]);
        }
    }
    else
    {
        #pragma omp parallel for schedule(static)
        for (int ik = 0; ik < kpoint_num; ik++)
        {
            diagonalize_GeneralizedSelfAdjointMatrix_range_1k(H_k[ik], S_k[ik], lower_eigen_index, upper_eigen_index, eigenvectors[ik], eigenvalues[ik]);
        }
    }

}