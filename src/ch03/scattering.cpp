// Chapter 3: Quantum mechanical scattering
// Single particle, one dimension, incident from left.
// Piecewise constant barriers, 0 elsewhere.
// (-h^2/2m x'' + v)w(x) = E*w(x)

// input file
// N            (N > 0 number of barriers)
// v1, v2, ...  (N barrier values)
// a1, a2, ...  (N+1 barrier boundaries)
// k            (sqrt(2mE)/hbar)
// xmin         (left endpoint of solution)
// xmax         (right endpoint of solution)
// M            (number of points in solution)

#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "futil/futil.h"

using Array = Eigen::ArrayXd;
using Vector = Eigen::VectorXd;
using CArray = Eigen::ArrayXcd;
using CVector = Eigen::VectorXcd;
using CMatrix = Eigen::MatrixXcd;

typedef std::complex<double> Complex;
typedef Eigen::Triplet<Complex> Triplet;

constexpr Complex I(0.0, 1.0);

// EDIT THE NAMESPACE FOR PIECEWISE FUNCTION GENERATION -g
namespace UserDef
{
    int N{ 10 };
    double amin{ -8.0 };
    double amax{ 8.0 };

    Array userFunc(const Array& x)
    {
        constexpr double a{ 2.0 };
        constexpr double v0{ 3.0 };
        return v0 * (x.square() / (-2.0 * a * a)).exp();
    }
}
// END OF EDIT

namespace Ex3_9
{
    Array gaussian(const Array& x, double a, double v0)
    {
        return v0 * (x.square() / (-2.0 * a * a)).exp();
    }
}

std::tuple<CArray, CVector> scattering(
    int N, const Array& v, const Array& a, double k)
{
    // compute entries of the matrix
    CArray n{ (CArray::Ones(N) - v).sqrt() };
    Complex ik{ I * k };
    CArray ink{ ik * n };

    CArray einka01{ (ink * a.head(N)).exp() };
    CArray einka01m{ (-1.0 * ink * a.head(N)).exp() };
    CArray einka11{ (ink * a.tail(N)).exp() };
    CArray einka11m{ (-1.0 * ink * a.tail(N)).exp() };

    CArray inkeinka01{ ink * einka01 };
    CArray inkeinka01m{ ink * einka01m };
    CArray inkeinka11{ ink * einka11 };
    CArray inkeinka11m{ ink * einka11m };

    Complex eika0{ std::exp(ik * a(0)) };
    Complex eika0m{ std::exp(-1.0 * ik * a(0)) };
    Complex eikaN{ std::exp(ik * a(N)) };

    Complex ikeika0{ ik * eika0 };
    Complex ikeika0m{ ik * eika0m };
    Complex ikeikaN{ ik * eikaN };

    // construct matrix
    // assume incoming wave has amplitude 1
    std::vector<Triplet> entries;
    entries.reserve(8*N+4);

    entries.push_back(Triplet(0, 0, eika0m));
    entries.push_back(Triplet(N, 2*N+1, -eikaN));
    entries.push_back(Triplet(N+1, 0, -ikeika0m));
    entries.push_back(Triplet(2*N+1, 2*N+1, -ikeikaN));

    for (int i{ 0 }; i < N; ++i)
    {
        entries.push_back(Triplet(i, 2*i+1, -einka01(i)));
        entries.push_back(Triplet(i, 2*i+2, -einka01m(i)));
        entries.push_back(Triplet(i+1, 2*i+1, einka11(i)));
        entries.push_back(Triplet(i+1, 2*i+2, einka11m(i)));

        entries.push_back(Triplet(i+N+1, 2*i+1, -inkeinka01(i)));
        entries.push_back(Triplet(i+N+1, 2*i+2, inkeinka01m(i)));
        entries.push_back(Triplet(i+N+2, 2*i+1, inkeinka11(i)));
        entries.push_back(Triplet(i+N+2, 2*i+2, -inkeinka11m(i)));
    }

    Eigen::SparseMatrix<Complex> A(2*N+2, 2*N+2);
    A.setFromTriplets(entries.begin(), entries.end());

    // construct rhs
    CVector b{ CVector::Zero(2*N+2) };
    b(0) = -eika0;
    b(N+1) = -ikeika0;

    // solve system
    Eigen::SparseLU<Eigen::SparseMatrix<Complex>> linsolve;
    linsolve.compute(A);
    CVector x{ linsolve.solve(b) };

    return {n, x};
}

std::tuple<Vector, CVector> constructSolution(
    const CArray& n, const CVector& coef, const Array& a,
    double k, double xmin, double xmax, int M)
{
    Vector x{ Vector::LinSpaced(M, xmin, xmax) };
    CVector psi(M);

    // first region
    int i{ 0 };
    for (; x(i) < a(0); ++i)
    {
        psi(i) = std::exp(I * k * x(i)) + coef(0) * std::exp(-I * k * x(i));
    }

    // regions in between
    for (int j{ 1 }; j < a.size(); ++j)
    {
        for (; x(i) < a(j); ++i)
        {
            psi(i) = coef(2*j-1) * std::exp(I * k * n(j-1) * x(i))
                + coef(2*j) * std::exp(-I * k * n(j-1) * x(i));
        }
    }

    // final region
    for (; i < M; ++i)
    {
        psi(i) = coef.tail(1)(0) * std::exp(I * k * x(i));
    }

    return {x, psi};
}

std::tuple<Vector, double, double> computeMisc(
    const CVector& psi, const CVector& coef, bool show)
{
    // probability density
    Vector prob{ psi.cwiseAbs2() };

    // reflection/transmission coefficients
    double refl{ std::norm(coef(0)) };
    double tran{ std::norm(coef.tail(1)(0)) };

    if (show)
    {
        std::cout << "Reflection coefficient: " << refl << '\n';
        std::cout << "Transmission coefficient: " << tran << '\n';
    }

    return {prob, refl, tran};
}

// converts a function to piecewise step function input file
void funcToStep(const std::string& outfile, std::function<Array(Array)> func,
    int N, double amin, double amax)
{
    Array a{ Array::LinSpaced(N+1, amin, amax) };
    Array v{ func(a.head(N)) };

    // write to file
    Eigen::IOFormat fmt(12, Eigen::DontAlignCols, " ", ",");
    std::ofstream output{ outfile };
    output << N << '\n' << v.format(fmt) << '\n' << a.format(fmt);
}

// Exercise 3.9
void exercise3_9(const std::string& outfile, double a, double v0, double k)
{
    int tot{ 14 };
    Vector nvec(tot);
    Vector tran(tot);
    for (int i{ 0 }; i < tot; ++i)
    {
        int N{ static_cast<int>(pow(2, i+1)) };
        nvec(i) = static_cast<double>(i+1);

        Array ai{ Array::LinSpaced(N+1, -4.0*a, 4.0*a) };
        Array vi{ Ex3_9::gaussian(ai.head(N), a, v0) };

        auto [n, coef] = scattering(N, vi, ai, k);
        tran(i) = std::norm(coef.tail(1)(0));
    }

    writeFile(outfile, 12, nvec, tran);
}

int main(int argc, const char* argv[])
{
    // default file names
    std::string infile{ "input.txt" };
    std::string outfile{ "output.txt" };
    std::string opflag{ "compute" };

    // read file names from command line
    for (int i{ 1 }; i < argc; i += 2)
    {
        if (!argv[i+1])
        {
            std::cerr << "Missing argument." << '\n';
            return 1;
        }

        std::string fileflag{ argv[i] };
        if (fileflag == "-i")
            infile = argv[i+1];
        else if (fileflag == "-o")
            outfile = argv[i+1];
        else if (fileflag == "-g")
        {
            outfile = argv[i+1];
            opflag = "generate";
        }
        else if (fileflag == "-3.9")
        {
            outfile = argv[i+1];
            opflag = "3.9";
        }
        else
        {
            std::cerr << fileflag << " is not valid." << '\n';
            return 1;
        }
    }

    // operation to perform
    if (opflag == "compute")
    {
        auto [N, v, a, k, xmin, xmax, M] = readFile
            <int, Array, Array, double, double, double, int>(infile);
        auto [n, coef] = scattering(N, v, a, k);
        auto [x, psi] = constructSolution(n, coef, a, k, xmin, xmax, M);
        auto [prob, refl, tran] = computeMisc(psi, coef, true);
        writeFile(outfile, 12, x, psi.real(), psi.imag(), prob);
    }
    else if (opflag == "generate")
    {
        funcToStep(outfile,
            UserDef::userFunc, UserDef::N, UserDef::amin, UserDef::amax);
    }
    else
    {
        exercise3_9(outfile, 0.5, 100.0, 0.1);
    }

    return 0;
}
