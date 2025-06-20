// Chapter 3: Coupled oscillators
// n masses connected by springs, n > 1
// mx'' = -kx

// input file
// n            (number of masses)
// m1, m2, ...  (n masses)
// x0, ...      (n initial displacements)
// v0, ...      (n initial velocities)
// k1, k2, ...  (n-1 spring constants)
// t0           (initial time)
// tN           (final time)
// N            (number of time steps)

// output file
// t, x1, ...   (times, displacement)

#include <fstream>
#include <iostream>
#include <tuple>

#include <Eigen/Dense>

#include "futil/futil.h"

using Array = Eigen::ArrayXd;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using CArray = Eigen::ArrayXcd;

constexpr std::complex<double> I(0.0, 1.0);
constexpr double sqrtzerotol{ 1e-7 };

std::tuple<Matrix, CArray, Array> coupledOscillators(
    int n, const Array& m, const Vector& x0, const Vector& v0, const Array& k)
{
    // construct K
    Matrix K(n, n);
    if (n > 2)
    {
        Vector kd(n);
        kd(0) = k(0);
        kd.segment(1, n-2) = k.head(n-2) + k.tail(n-2);
        kd(n-1) = k(n-2);

        K.diagonal() << kd;
        for (int i{ 0 }; i < n-1; ++i)
        {
            K(i, i+1) = -k(i);
            K(i+1, i) = -k(i);
        }
    }
    else
        K << k(0), -k(0), -k(0), k(0);

    // construct M^(-1/2)
    Vector m12{ m.sqrt() };
    Vector m12inv{ m12.cwiseInverse() };
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M12inv{ m12inv };

    // solve eigenvalue problem
    Matrix A{ M12inv * K * M12inv };
    Eigen::SelfAdjointEigenSolver<Matrix> eigsol(A);
    Array eigvals{ eigsol.eigenvalues().cwiseSqrt() };
    Matrix eigvecs{ eigsol.eigenvectors() };

    // compute coefficients
    Matrix coef_vect(n, n);
    CArray coef_scal(n);
    for (int i{ 0 }; i < n; ++i)
    {
        coef_vect.col(i) = eigvecs.col(i).cwiseProduct(m12inv);

        // bT*M
        Eigen::RowVectorXd btm{ eigvecs.col(i).cwiseProduct(m12).transpose() };

        // check for zero eigenvalues
        if (std::abs(eigvals(i)) < sqrtzerotol || std::isnan(eigvals(i)))
        {
            coef_scal(i) = btm * x0;
            eigvals(i) = 0.0;
        }
        else
            coef_scal(i) = (btm * (x0 + I / eigvals[i] * v0))(0);
    }

    return {coef_vect, coef_scal, eigvals};
}

std::tuple<Vector, Matrix> constructSolution(
    const Matrix& a, const CArray& c, const Array& w,
    double t0, double tN, int N)
{
    // time grid
    Vector t{ Vector::LinSpaced(N+1, t0, tN) };

    // solution
    Matrix x(N+1, c.size());
    for (int i{ 0 }; i < N+1; ++i)
    {
        Eigen::VectorXcd ceiwt{ c * (-1.0 * I * t[i] * w).exp() };
        x.row(i) = (a * ceiwt).real();
    }

    return {t, x};
}

int main(int argc, const char* argv[])
{
    // default file names
    std::string infile{ "input.txt" };
    std::string outfile{ "output.txt" };

    // read file names from command line
    for (int i{ 1 }; i < argc; i += 2)
    {
        if (!argv[i+1])
        {
            std::cerr << "Missing argument." << '\n';
            return 1;
        }

        std::string flag{ argv[i] };
        if (flag == "-i")
            infile = argv[i+1];
        else if (flag == "-o")
            outfile = argv[i+1];
        else
        {
            std::cerr << flag << " is not valid." << '\n';
            return 1;
        }
    }

    auto [n, m, x0, v0, k, t0, tN, N] = readFile
        <int, Array, Vector, Vector, Array, double, double, int>(infile);
    auto [a, c, w] = coupledOscillators(n, m, x0, v0, k);
    auto [t, x] = constructSolution(a, c, w, t0, tN, N);
    writeFile(outfile, 12, t, x);

    return 0;
}
