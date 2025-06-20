#include "futil/futil.h"
#include "mutil/interp.h"

#include <iostream>
#include <iomanip>

using Array = Eigen::ArrayXd;

// Exercise 4.2
void exercise4_2(int n, const std::string& outf)
{
    // data
    const Eigen::Array<double, 1, 7> xd{ -1.0, -0.5*sqrt(3.0), -0.5, 0.0, 0.5, 0.5*sqrt(3.0), 1.0 };
    const Eigen::Array<double, 1, 7> yd{ 0.0, 0.5, 0.5*sqrt(3.0), 1.0, 0.5*sqrt(3.0), 0.5, 0.0 };

    // interpolate
    const auto [x, y_lagrange] = lagrange(xd, yd, n);
    const Array y_spline1{ std::get<1>(cubic_spline1(xd, yd, n)) };
    const Array y_spline2{ std::get<1>(cubic_spline2(xd, yd, n)) }; 
    const Array y_exact{ (1.0 - x*x).sqrt() };

    writeFile(outf, 12, x, y_exact, y_lagrange, y_spline1, y_spline2);
}

// Exercise 4.7
void exercise4_7(int n, const std::string& inf, const std::string& outf)
{
    const auto [xd, yd] = readFile<Array, Array>(inf);
    const auto [x, y_spline1] = cubic_spline1(xd, yd, n);
    writeFile(outf, 12, x, y_spline1);
}

// Exercise 4.8
void exercise4_8(int n, const std::string& inf, const std::string& outf)
{
    const auto [xblu, yblu, xrod, yrod, xgrn, ygrn, xred, yred] = readFile
        <Array, Array, Array, Array, Array, Array, Array, Array>(inf);

    const auto [xintblu, yintblu] = cubic_spline1(xblu, yblu, n);
    const auto [xintrod, yintrod] = cubic_spline1(xrod, yrod, n);
    const auto [xintgrn, yintgrn] = cubic_spline1(xgrn, ygrn, n);
    const auto [xintred, yintred] = cubic_spline1(xred, yred, n);

    writeFile(outf, 12, xintblu, yintblu, xintrod, yintrod, xintgrn, yintgrn, xintred, yintred);
}

// Exercise 4.9
void exercise4_9(const std::string& inf, const std::string& outf)
{
    // initialize
    const auto [xblu, yblu, _1, _2, xgrn, ygrn, xred, yred] = readFile
        <Array, Array, Array, Array, Array, Array, Array, Array>(inf);
    Array yintgrn{ Array::Zero(28001) };
    Array yintblu{ Array::Zero(28001) };

    // interpolate
    const auto [xint, yintred] = cubic_spline1(xred, yred, 28001);
    yintgrn.head(25001) = std::get<1>(cubic_spline1(xgrn, ygrn, 25001));
    yintblu.head(13001) = std::get<1>(cubic_spline1(xblu.tail(14), yblu.tail(14), 13001));

    // concatenate, cap at 100, and normalize
    Eigen::ArrayXXd rgb
    { 
        (Eigen::ArrayXXd(28001, 3) << yintred, yintgrn, yintblu).finished()
        .min(100.0) * 2.55
    };

    // create image
    int n{ 101 };
    int nrows{ n*28001 };
    Eigen::ArrayXXd im(nrows, 5);
    Array y{ Array::LinSpaced(n, 0.0, 1.0) };
    for (int i{ 0 }; i < nrows; ++i)
    {
        im.row(i) << xint(i / n), y(i % n), rgb(i / n, 0), rgb(i / n, 1), rgb(i / n, 2);
    }

    writeFile(outf, 12, im);
}

// Exercise 4.10
void exercise4_10()
{
    // exact answer
    double x{ 0.99 };
    double f_exact{ 1.0 / ((1.0 + x) * (2.0 + x)) };
    std::cout << std::setprecision(9) << "f(0.99):\t" << f_exact << '\n';

    // initialize (i = 0)
    int i{ 0 };
    double x_pow{ 1.0 };
    double half_pow{ 0.5 };
    double alt_term{ 1.0 };
    double f_approx{ 0.5 };

    // Taylor aproximation
    while (abs(f_exact - f_approx) > 1e-7)
    {
        x_pow *= 0.99;
        half_pow *= 0.5;
        alt_term = -alt_term;

        f_approx += (1.0 - half_pow) * x_pow * alt_term;
        ++i;
    }
    std::cout << std::setprecision(9) << "i = " << i << ": \t" << f_approx << '\n';
}

// Exercise 4.11
void exercise4_11()
{
    // get input
    std::cout << "Enter a double with absolute value less than 1: ";
    double x{};
    std::cin >> x;

    std::cout << "Enter an integer greater than 1: ";
    int n{};
    std::cin >> n;

    double n_partial{ (1.0 - pow(-x, n+1)) / (1.0 + x) };
    double nplus_partial{ (1.0 - pow(-x, n+2)) / (1.0 + x) };
    double nminus_partial{ (1.0 - pow(-x, n)) / (1.0 + x) };
    double shanks{ (nplus_partial * nminus_partial - n_partial * n_partial)
        / (nplus_partial + nminus_partial - 2.0 * n_partial) };
    
    double series_sum{ 1.0 / (1.0 + x) };

    std::cout << "The sum of the series: \t" << series_sum << '\n';
    std::cout << "Shanks transformation: \t" << shanks << '\n';
}

// Exercise 4.12
void exercise4_12(const int ord, const int n)
{
    // first partial sums
    auto sums = [](int n)
    {
        std::vector<double> partials(n);
        double alt_term{ 1.0 };
        double psum{};

        for (int i{ 0 }; i < n; ++i)
        {
            alt_term = -alt_term;
            psum += alt_term / sqrt(i+1);
            partials[i] = psum;
        }

        return partials;
    };

    std::cout << std::setprecision(9) << shanks(sums, ord, n) << '\n';
}

// Exercise 4.13
void exercise4_13(const std::string& outf)
{
    int n{ 30 };
    Array pi_vals(n);

    for (int i{ 1 }; i <= n; ++i)
    {
        Array a_vals(i);
        a_vals(0) = 1;
        a_vals.tail(i-1) = Array::LinSpaced(i-1, 1, i-1).square();

        Array b_vals{ Array::LinSpaced(i, 1, 2*i-1) };

        pi_vals(i-1) = 4.0 * cf_to_pade1(a_vals, b_vals);
    }

    writeFile(outf, 12, Array::LinSpaced(n, 1, n), (M_PI - pi_vals).abs());
}

// Exercise 4.14
void exercise4_14()
{   
    const int n{ 8 };
    const double alpha{ 1.0 };
    const double beta{ 1.0 };

    // matrix values
    const double beta2{ beta * beta };
    const Array m{ Array::LinSpaced(n, 0, n-1) };

    const Array diag{ m*(m + 2.0) + 2.0 * beta2 };
    const Array offdiag{ -beta2 * Array::Ones(n-1) };

    // convert linear system to continued fraction
    // f(x) = 0: the roots of the function gives the eigenvalues
    const Array a_vals{ - offdiag.square() };
    auto f = [diag, a_vals](double x)
    {
        Array b_vals{ x - diag.tail(n-1) };
        auto [num, denom] = cf_to_pade2(a_vals, b_vals);
        return (x - diag(0)) * denom + num;
    };

    // find zeros (on range given by Gershgorin circles)
    std::vector zeros{ find_zeros(f, diag(0)-beta2, diag(n-1)+beta2) };

    for (double val : zeros)
    {
        std::cout << std::setprecision(12) << val << '\n';
    }
}

// Exercise 4.15
void exercise4_15()
{
    
}



int main(const int argc, const char* argv[])
{
    // exercise4_2(1001, "out.txt");
    // exercise4_7(801, "input_4-7.txt", "out.txt");
    // exercise4_8(1001, "input_4-8.txt", "out.txt");
    // exercise4_9("input_4-8.txt", "out.txt");
    // exercise4_10();
    // exercise4_11();
    // exercise4_12(4, 5);
    // exercise4_13("out.txt");
    // exercise4_14();
    // exercise4_15();




    // for (int i{ 1 }; i < argc; ++i)
    // {

    // }
    // switch(*argv[1])
    // {
    //     case '2': exercise4_2(201, "out.txt");
    //     default: break;
    // }


    return 0;
}
