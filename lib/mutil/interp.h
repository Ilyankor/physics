#pragma once

#include <tuple>
#include <vector>

#include <Eigen/Dense>

using Array = Eigen::ArrayXd;

// numerical computation of Lagrange interpolation polynomial
// using Neville's algorithm
std::tuple<Array, Array> lagrange(const Array& xd, const Array& yd, const int N);

// cubic spline interpolation
// natural parameters
// first derivative based
std::tuple<Array, Array> cubic_spline1(const Array& xd, const Array& yd, const int N);

// cubic spline interpolation
// natural parameters
// second derivative based
std::tuple<Array, Array> cubic_spline2(const Array& xd, const Array& yd, const int N);

// shanks transformation to compute sum of series
double shanks(std::function<std::vector<double>(int)> f, const int ord, const int n);

// convert continued fraction to pade approximation
double cf_to_pade1(const Array& a_vals, const Array& b_vals);

// returns numerator and denominator
std::tuple<double, double> cf_to_pade2(const Array& a_vals, const Array& b_vals);

// bisection rootfinding algorithm
double bisection(std::function<double(double)> f, double a, double b, 
    const double tol = 1e-12, const int max_iter = 1000, const double zerotol = 1e-14);

// find zeros
std::vector<double> find_zeros(std::function<double(double)> f,
    const double x_start, const double x_end, const double n = 1000,
    const double tol = 1e-12, const int max_iter = 1000, const double zerotol = 1e-14);

// pade approximation of a power series
double series_to_pade(const Array& a, double x);
