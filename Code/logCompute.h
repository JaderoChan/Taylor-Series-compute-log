#pragma once

#include <cmath>

#include <iostream>
#include <exception>

namespace MyLog
{

/*
$$
\because ln(2^\frac{1}{16})=\frac{1}{16}ln(2)
\\
\therefore ln(2)=16ln(2^\frac{1}{16})
\\
\ast \quad 2-1=1, \ \sqrt[16]{2}-1\approx 0.0442737 \quad \ast
$$
*/
// @brief Compute the ln(2) by Taylor series.
// @param series The series larger, result more precise, but no obvious effects if the series over 7 places.
inline double ln2ByTaylorSeries(int series = 1'000'000) {
    double sum = 0;
    double x = std::pow(2, 0.0625) - 1;
    for (int i = 1; i <= series; ++i) {
        if (i % 2 == 0)
            sum -= std::pow(x, i) / i;
        else
            sum += std::pow(x, i) / i;
    }
    return 16 * sum;
}

/*
$$
\because ln(3^\frac{1}{16})=\frac{1}{16}ln(3)
\\
\therefore ln(3)=16ln(3^\frac{1}{16})
\\
\ast \quad 3-1=2, \ \sqrt[16]{3}-1\approx 0.0710754 \quad \ast
$$
*/
// @brief Compute the ln(2) by Taylor series.
// @param series The series larger, result more precise, but no obvious effects if the series over 7 places.
inline double ln3ByTaylorSeries(int series = 1'000'000) {
    double sum = 0;
    double x = std::pow(3, 0.0625) - 1;
    for (int i = 1; i <= series; ++i) {
        if (i % 2 == 0)
            sum -= std::pow(x, i) / i;
        else
            sum += std::pow(x, i) / i;
    }
    return 16 * sum;
}

const double Ln2 = ln2ByTaylorSeries();
const double Ln3 = ln3ByTaylorSeries();

// @param series The series larger, result more precise.
// @note The m in range[1,2), 2m/3 in range[2/3, 4/3). Let "x" moer close 1 for higher precision.
inline double ln(double x, int series = 10) {
    if (x <= 0)      // The log() domain of definition is (0, +infty).
        throw std::domain_error("The x must is positive.");
    if (x == 1)    // The ln(1) = 0.
        return 0;
    // Cast the double to unsigned long long for bit operation.
    unsigned long long _x = *reinterpret_cast<unsigned long long *>(&x);
    // Get the float point num's exp.
    int j = _x >> 52;
    // Get the float point num's mantissa. The (unsigned long long(0x3ff0) << 48) is a double it = 1.0.
    _x = unsigned long long(0x3ff0) << 48 ^ (_x << 12 >> 12);
    // Cast the bit operated unsigned long long to double.
    double m = *reinterpret_cast<double *>(&_x);
    // Let "x" more close 1 for higher precision.
    m = 2 * m / 3;
    // The Taylor series expanded summation.
    double sum = 0;
    for (int i = 1; i <= series; ++i) {
        if (i % 2 == 0)
            sum -= std::pow((m - 1), i) / i;
        else
            sum += std::pow((m - 1), i) / i;
    }
    return sum + Ln3 - Ln2 + (j - 1023) * Ln2;
}

// @param series The series larger, result more precise
// @note The m in range[1,2), √2m/2 in range[√2/2, √2). Let x moer close 1 for get higher precision.
inline double ln_(double x, int series = 10) {
    if (x <= 0)
        throw std::domain_error("The x must is positive.");
    if (x == 1)
        return 0;
    unsigned long long _x = *reinterpret_cast<unsigned long long *>(&x);
    int j = _x >> 52;
    _x = unsigned long long(0x3ff0) << 48 ^ (_x << 12 >> 12);
    double m = *reinterpret_cast<double *>(&_x);
    m = std::pow(2, 0.5) * m / 2;
    double sum = 0;
    for (int i = 1; i <= series; ++i) {
        if (i % 2 == 0)
            sum -= std::pow((m - 1), i) / i;
        else
            sum += std::pow((m - 1), i) / i;
    }
    return sum + Ln2 / 2 + (j - 1023) * Ln2;
}

inline double log(double base, double x) {
    if (base == 1)
        throw std::domain_error("The base num can't is 1.");
    return ln(x) / ln(base);
}

}

namespace StdLog
{

inline double ln(double x) {
    return std::log(x);
}

inline double log(double base, double x) {
    return ln(x) / ln(base);
}

}

namespace UnitTest
{

inline void logPrecisionCompare(int testCount = 20) {
    std::cout.precision(12);

    double x = 0.0;
    for (int i = 0; i < testCount; ++i) {
        x += 0.2;
        double a1 = MyLog::ln(x);
        double a2 = StdLog::ln(x);

        std::cout << "X: " << x << ";   \t";
        std::cout << "My: " << a1 << ";   \t";
        std::cout << "STD: " << a2 << ";   \t";
        std::cout << "Deviation: " << a2 - a1 << "\n";
    }
}

inline void logPerformanceCompare(int testCount = 100'000) {
    time_t start = clock();

    double x = 0.25;
    for (int i = 0; i < testCount; ++i)
        StdLog::ln(x++);

    std::cout << "[STD][" << testCount << "]: ";
    std::cout << clock() - start << "msec" << '\n';
    start = clock();

    x = 0.25;
    for (int i = 0; i < testCount; ++i)
        MyLog::ln(x++);

    std::cout << "[MY][" << testCount << "]: ";
    std::cout << clock() - start << "msec" << '\n';
}

}
