#include "futil.h"
#include <sstream>
#include <vector>

Array strToArr(const std::string& str)
{
    std::vector<double> values;
    std::stringstream ss{ str };

    double val{};
    char sep{};
    while (ss >> val)
    {
        values.push_back(val);
        ss >> sep;
    }

    return Eigen::Map<Array>(values.data(), values.size());
}
