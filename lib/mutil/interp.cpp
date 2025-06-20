#include "interp.h"

#include <iostream>

std::tuple<Array, Array> lagrange(const Array& xd, const Array& yd, const int N)
{
    // create solution points
    const int M{ static_cast<int>(xd.size()) };
    Array x{ Array::LinSpaced(N, xd(0), xd(M-1)) };
    
    // initialize
    std::vector<Array> p(M);
    for (int i{ 0 }; i < M; ++i) { p[i] = (yd(i) * Array::Ones(N)); }

    // Neville algorithm
    for (int i{ 1 }; i < M; ++i)
    {
        for (int j{ 0 }; j < M-i; ++j)
        {
            p[j] = ((x - xd(j)) * p[j+1] - (x - xd(j+i)) * p[j])
                / (xd(j+i) - xd(j));
        }
    }

    return std::tuple{ x, p[0] };
}

std::tuple<Array, Array> cubic_spline1(const Array& xd, const Array& yd, const int N)
{
    // create arrays
    const int M{ static_cast<int>(xd.size()) - 1 };
    const Array dx{ xd.tail(M) - xd.head(M) };
    const Array dy{ yd.tail(M) - yd.head(M) };

    const Array invdx{ dx.inverse() };
    const Array m{ dy * invdx };
    const Array mdx{ 3.0 * m * invdx };

    // create matrix
    Eigen::VectorXd diag(M+1);
    diag(0) = 2.0 * invdx(0);
    diag(M) = 2.0 * invdx(M-1);
    diag.segment(1, M-1) = 2.0 * (invdx.head(M-1) + invdx.tail(M-1));

    Eigen::MatrixXd A{ diag.asDiagonal() };
    for (int i{ 0 }; i < M; ++i)
    {
        A(i, i+1) = invdx(i);
        A(i+1, i) = invdx(i);
    }

    // construct rhs
    Eigen::VectorXd rhs(M+1);
    rhs(0) = mdx(0);
    rhs(M) = mdx(M-1);
    rhs.segment(1, M-1) = mdx.head(M-1) + mdx.tail(M-1);

    // linear system
    Eigen::VectorXd q{ A.colPivHouseholderQr().solve(rhs) };

    // construct solution
    int i{ 0 };
    Array x{ Array::LinSpaced(N, xd(0), xd(M)) };
    Array sol(N);
    for (int j{ 0 }; j < M; ++j)
    {
        for (; x(i) < xd(j+1); ++i)
        {
            double t1{ (x(i) - xd(j)) / dx(j) };
            double t2{ 1.0 - t1 };
            sol(i) = t2 * yd(j) + t1 * yd(j+1) + t1 * t2 * dx(j)
                * ((m(j) - q(j+1)) * t1 + (q(j) - m(j)) * t2);
        }
    }
    sol.tail(1) = yd.tail(1);

    return std::tuple{ x, sol };
}

std::tuple<Array, Array> cubic_spline2(const Array& xd, const Array& yd, const int N)
{
    // create arrays
    const int M{ static_cast<int>(xd.size()) - 1 };

    const Array dx{ xd.tail(M) - xd.head(M) };
    const Array dy{ yd.tail(M) - yd.head(M) };
    const Array m{ dy / dx };

    // create matrix
    Eigen::VectorXd diag{ (dx.segment(0, M-1) +  dx.segment(1, M-1)) / 3.0 };
    Eigen::MatrixXd A{ diag.asDiagonal() };
    for (int i{ 0 }; i < M-2; ++i)
    {
        A(i, i+1) = dx(i+1) / 6.0;
        A(i+1, i) = dx(i+1) / 6.0;
    }

    // construct rhs
    Eigen::VectorXd rhs{ m.tail(M-1) - m.head(M-1) };

    // linear system
    Array z{ Array::Zero(M+1) };
    z.segment(1, M-1) = A.colPivHouseholderQr().solve(rhs);

    // construct solution
    int i{ 0 };
    Array x{ Array::LinSpaced(N, xd(0), xd(M)) };
    Array sol(N);
    for (int j{ 0 }; j < M; ++j)
    {
        for (; x(i) < xd(j+1); ++i)
        {
            double t1{ x(i) - xd(j) };
            double t2{ x(i) - xd(j+1) };
            sol(i) = (z(j+1) * pow(t1, 3) - z(j) * pow(t2, 3)
                + (6.0 * yd(j+1) - z(j+1) * pow(dx(j), 2)) * t1
                - (6.0 * yd(j) - z(j) * pow(dx(j), 2)) * t2)
                / (6.0 * dx(j));
        }
    }
    sol.tail(1) = yd.tail(1);

    return std::tuple{ x, sol };
}

double shanks(std::function<std::vector<double>(int)> f, const int ord, const int n)
{
    // initialize
    int nterms{ n + 2*(ord - 1) };
    std::vector<double> partials{ f(nterms) };

    // next order partial sums
    for (int i{ 0 }; i < ord; ++i)
    {
        for (int j{ 1 }; j < nterms-1-2*i; ++j)
        {
            partials[j-1] = (partials[j+1] * partials[j-1] - partials[j] * partials[j])
                / (partials[j+1] + partials[j-1] - 2.0 * partials[j]);
        }
    }

    return partials[n-1];
}

double cf_to_pade1(const Array& a_vals, const Array& b_vals)
{
    double a_prev{ 1.0 };
    double a_cur{ 0.0 };

    double b_prev{ 0.0 };
    double b_cur{ 1.0 };

    for (int i{ 0 }; i < a_vals.size(); ++i)
    {
        double a_new{ b_vals(i) * a_cur + a_vals(i) * a_prev };
        double b_new{ b_vals(i) * b_cur + a_vals(i) * b_prev };

        a_prev = a_cur;
        a_cur = a_new;

        b_prev = b_cur;
        b_cur = b_new;
    }

    return a_cur / b_cur;
}

std::tuple<double, double> cf_to_pade2(const Array& a_vals, const Array& b_vals)
{
    double a_prev{ 1.0 };
    double a_cur{ 0.0 };

    double b_prev{ 0.0 };
    double b_cur{ 1.0 };

    for (int i{ 0 }; i < a_vals.size(); ++i)
    {
        double a_new{ b_vals(i) * a_cur + a_vals(i) * a_prev };
        double b_new{ b_vals(i) * b_cur + a_vals(i) * b_prev };

        a_prev = a_cur;
        a_cur = a_new;

        b_prev = b_cur;
        b_cur = b_new;
    }

    return std::tuple(a_cur, b_cur);
}

double bisection(std::function<double(double)> f,
    double a, double b,
    const double tol, const int max_iter, const double zerotol)
{
    double p{};
    double fp{};
    double len{ b - a };

    for (int i{ 0 }; i < max_iter; ++i)
    {
        p = 0.5 * (a + b);
        fp = f(p);

        // check if 0 is reached
        if (abs(fp) < zerotol) { return p; }
        
        // check if tolerance is reach
        len *= 0.5;
        if (len < tol) { return p; }

        // next iteration
        if (f(a) * fp < 0.0) { b = p; }
        else { a = p; }
    }

    return p;
}

std::vector<double> find_zeros(std::function<double(double)> f,
    const double x_start, const double x_end, const double n,
    const double tol, const int max_iter, const double zerotol)
{
    // sample points
    const Array x_sample{ Array::LinSpaced(n+1, x_start, x_end) };
    const Array a_vals{ x_sample.head(n) };
    const Array b_vals{ x_sample.tail(n) };

    // bisection method
    std::vector<double> zeros{};
    for (int i{ 0 }; i < n; ++i)
    {
        double a{ a_vals(i) };
        double b{ b_vals(i) };

        double fa{ f(a) };
        double fb{ f(b) };

        if (abs(fa) < zerotol) { zeros.push_back(a); }
        else if (abs(fb) < zerotol) { zeros.push_back(b); }
        else if (fa * fb < 0) { zeros.push_back(
            bisection(f, a_vals(i), b_vals(i), tol, max_iter, zerotol)); }
    }

    return zeros;
}

std::vector<double> quotient_difference(const Array& a)
{
    // initialize
    const int n = a.size() - 1;
    Array e{ Array::Zero(n) };
    Array q{ a.tail(n) / a.head(n) };

    std::vector<double> sol;
    sol.push_back(a(0));
    sol.push_back(-q(0));

    // quotient difference algorithm
    for (int i{ 1 }; i <= n/2; ++i)
    {
        e = q.tail(n-i) - q.head(n-i) + e.tail(n-i).eval();
        sol.push_back(-e(0));

        q = e.tail(n-i-1) / e.head(n-i-1) * q.segment(1, n-i-1).eval();
        sol.push_back(-q(0));
    }

    sol.pop_back();
    return sol;
}

double series_to_pade(const Array& a, double x)
{
    std::vector<double> c{ quotient_difference(a) };
    
    double pade{ 1.0 };
    for (int i = c.size()-1; i > 0; --i)
        pade = 1.0 + c[i]*x / pade;

    return c[0] / pade;
}
