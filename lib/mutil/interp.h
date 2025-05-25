#pragma once

#include <vector>

#include <Eigen/Core>

using Array = Eigen::ArrayXd;

// numerical computation of Lagrange interpolation polynomial
// using Neville's algorithm
Array lagrange(const Array& x,
    const std::vector<double>& xd, const std::vector<double>& yd);

double lagrange(double x,
    const std::vector<double>& xd, const std::vector<double>& yd);
