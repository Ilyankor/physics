#include "interp.h"

Array lagrange(const Array& x,
    const std::vector<double>& xd, const std::vector<double>& yd)
{
    int xlen = x.size();
    std::vector<Array> p{};
    for (auto& y : yd) { p.push_back(y * Array::Ones(xlen)); }

    int n = xd.size();
    for (int i{ 1 }; i < n; ++i)
    {
        for (int j{ 0 }; j < n-i; ++j)
        {
            p[j] = ((x - xd[j]) * p[j+1] - (x - xd[j+i]) * p[j])
                / (xd[j+i] - xd[j]);
        }
    }

    return p[0];
}

double lagrange(double x,
    const std::vector<double>& xd, const std::vector<double>& yd)
{
    std::vector<double> p{ yd };

    int n = xd.size();
    for (int i{ 1 }; i < n; ++i)
    {
        for (int j{ 0 }; j < n-i; ++j)
        {
            p[j] = ((x - xd[j]) * p[j+1] - (x - xd[j+i]) * p[j])
                / (xd[j+i] - xd[j]);
        }
    }

    return p[0];
}
