#include "util.h"
#include <sstream>

Eigen::ArrayXd strToArr(const std::string& str, int n)
{
    Eigen::ArrayXd arr(n);
    std::stringstream ss{ str };

    double val{};
    char sep{};

    for (int i{ 0 }; i < n; ++i)
    {
        ss >> val >> sep;
        arr(i) = val;
    }

    return arr;
}
