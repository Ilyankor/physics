// Chapter 3: Coupled oscillators
// n masses connected by springs, n > 1
// mx'' = -kx

// input file
// n            (number of masses)
// m1, m2, ...  (masses)
// x0, ...      (initial displacements)
// v0, ...      (initial velocities)
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

using Vector = Eigen::VectorXd;
using CVector = Eigen::VectorXcd;
using Matrix = Eigen::MatrixXd;

constexpr std::complex<double> I(0.0, 1.0);

Vector strToVec(std::string& str, int n)
{
    Vector vec(n);
    std::stringstream ss{ str };

    double val {};
    char comma {};

    for (int i{ 0 }; i < n; ++i)
    {
        ss >> val >> comma;
        vec[i] = val;
    }

    return vec;
}

std::tuple<int, Vector, Vector, Vector, Vector, double, double, int> readFile
    (char* file_name)
{
    std::ifstream input{ file_name };
    std::string line {};

    std::getline(input, line);
    int n{ std::stoi(line) };
    
    std::getline(input, line);
    Vector m{ strToVec(line, n) };

    std::getline(input, line);
    Vector x0{ strToVec(line, n) };

    std::getline(input, line);
    Vector v0{ strToVec(line, n) };

    std::getline(input, line);
    Vector k{ strToVec(line, n-1) };
    
    std::getline(input, line);
    double t0{ std::stod(line) };

    std::getline(input, line);
    double tN{ std::stod(line) };

    std::getline(input, line);
    int N{ std::stoi(line) };
    return {n, m, x0, v0, k, t0, tN, N};
}

void writeFile(char* file_name, Vector& t, Matrix& x)
{
    Matrix merged(x.rows(), x.cols() + 1);
    merged << t, x;

    Eigen::IOFormat fmt(12, 0, "\t", "\n");
    std::ofstream output{ file_name };
    output << merged.format(fmt);
}

std::tuple<Matrix, CVector, Vector> coupledOscillators
    (int n, Vector& m, Vector& x0, Vector& v0, Vector& k)
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
    {
        K << k(0), -k(0), -k(0), k(0);
    }

    // construct M^(-1/2)
    Vector m12{ m.cwiseSqrt() };
    Vector m12inv{ m12.cwiseInverse() };
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M12inv{ m12inv };

    // solve eigenvalue problem
    Matrix A{ M12inv * K * M12inv };
    Eigen::SelfAdjointEigenSolver<Matrix> eigsol(A);
    Vector eigvals{ eigsol.eigenvalues().cwiseSqrt() };
    Matrix eigvecs{ eigsol.eigenvectors() };

    // compute coefficients
    Matrix coef_vect(n, n);
    CVector coef_scal(n);
    for (int i{ 0 }; i < n; ++i)
    {
        coef_vect.col(i) = eigvecs.col(i).cwiseProduct(m12inv);

        // check zero eigenvalue
        if (std::abs(eigvals[i]) < 1e-7 || std::isnan(eigvals[i]))
        {
            coef_scal[i] = eigvecs.col(i).cwiseProduct(m12).transpose() * x0;
            eigvals[i] = 0.0;
        }
        else
        {
            coef_scal[i] = (eigvecs.col(i).cwiseProduct(m12).transpose() 
                * (x0 + I / eigvals[i] * v0)).value();
        }
    }

    return {coef_vect, coef_scal, eigvals};
}

std::tuple<Vector, Matrix> constructSolution
    (Matrix& a, CVector& c, Vector& w, double t0, double tN, int N)
{
    // time grid
    Vector t{ Vector::LinSpaced(N+1, t0, tN) };

    // solution
    Matrix x(N+1, c.rows());
    for (int i{ 0 }; i < N+1; ++i)
    {
        CVector ceiwt{ c.array() * (-1.0 * I * t[i] * w.array()).exp() };
        x.row(i) = (a * ceiwt).real();
    }

    return {t, x};
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cout << "Must have 2 arguments.\n";
        exit(1);
    }

    auto [n, m, x0, v0, k, t0, tN, N] = readFile(argv[1]);
    auto [a, c, w] = coupledOscillators(n, m, x0, v0, k);
    auto [t, x] = constructSolution(a, c, w, t0, tN, N);
    writeFile(argv[2], t, x);

    return 0;
}
