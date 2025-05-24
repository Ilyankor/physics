#pragma once

#include <fstream>
#include <string>
#include <Eigen/Core>

using Vector = Eigen::VectorXd;
using Array  = Eigen::ArrayXd;

// convert string to Eigen::ArrayXd
Array strToArr(const std::string& str);

// read a single line to various types
template <typename T>
inline T readLine(std::ifstream& input);

template <>
inline int readLine<int>(std::ifstream& input)
{
    std::string line{};
    std::getline(input, line);
    return std::stoi(line);
}

template <>
inline double readLine<double>(std::ifstream& input)
{
    std::string line{};
    std::getline(input, line);
    return std::stod(line);
}

template <>
inline Array readLine<Array>(std::ifstream& input)
{
    std::string line{};
    std::getline(input, line);
    return strToArr(line);
}

template <>
inline Vector readLine<Vector>(std::ifstream& input)
{
    std::string line{};
    std::getline(input, line);
    return strToArr(line).matrix();
}

// read input file
template <typename... Ts>
inline std::tuple<Ts...> readFile(const std::string& infile)
{
    std::ifstream stream{ infile };
    return std::make_tuple(readLine<Ts>(stream) ...);
}

// write data to file
template <typename... Ts>
inline void writeFile(std::string_view outfile, int precision, const Ts&... args)
{
    const int numrows = std::get<0>(std::make_tuple(args.rows() ...));
    const int totcols = (args.cols() + ...);
    Eigen::MatrixXd data(numrows, totcols);

    int col{ 0 };
    auto assign = [&](const auto& arr)
        {
            const int numcols = arr.cols();
            data.block(0, col, numrows, numcols) = arr;
            col += numcols;
        };
    (assign(args), ...);

    Eigen::IOFormat fmt(precision, 0, "\t", "\n");
    std::ofstream output{ outfile };
    output << data.format(fmt);
}
