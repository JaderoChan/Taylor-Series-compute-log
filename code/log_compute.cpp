#include <cmath>
#include <iostream>
#include <exception>

namespace my_log
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
/// @brief Compute the ln(2) by Taylor series.
/// @param series The series larger, result more precise, but no obvious effects if the series over 7 places.
double ln2ByTaylorSeries(int series = 1e6) {
    double sum = 0.0;
    double x = std::pow(2.0, 0.0625) - 1.0;
    bool isPositive = true;
    for (int i = 1; i <= series; ++i) {
        if (isPositive)
            sum += std::pow(x, i) / i;
        else
            sum -= std::pow(x, i) / i;
        isPositive = !isPositive;
    }
    return 16.0 * sum;
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
/// @brief Compute the ln(2) by Taylor series.
/// @param series The series larger, result more precise, but no obvious effects if the series over 7 places.
double ln3ByTaylorSeries(int series = 1e6) {
    double sum = 0.0;
    double x = std::pow(3.0, 0.0625) - 1.0;
    bool isPositive = true;
    for (int i = 1; i <= series; ++i) {
        if (isPositive)
            sum += std::pow(x, i) / i;
        else
            sum -= std::pow(x, i) / i;
        isPositive = !isPositive;
    }
    return 16.0 * sum;
}

const double LN2 = ln2ByTaylorSeries();
const double LN3 = ln3ByTaylorSeries();

/// @param series The series larger, result more precise.
/// @note The m in range[1,2), 2m/3 in range[2/3, 4/3). Let "x" moer close 1 for higher precision.
double lnVersion1(double x, int series = 10) {
    if (x <= 0.0)      // The log() domain of definition is (0, +infty).
        throw std::domain_error("The x must is positive.");
    if (x == 1.0)    // The ln(1) = 0.
        return 0.0;
    // Cast the double to unsigned long long for bit operation.
    unsigned long long _x = *reinterpret_cast<unsigned long long*>(&x);
    // Get the float point num's exp.
    int j = _x >> 52;
    // Get the float point num's mantissa. The (0x3ff0ULL << 48) is a double it = 1.0.
    _x = 0x3ff0ULL << 48 ^ (_x << 12 >> 12);
    // Cast the bit operated unsigned long long to double.
    double m = *reinterpret_cast<double*>(&_x);
    // Let "x" more close 1 for higher precision.
    m = 2.0 * m / 3.0;
    // The Taylor series expanded summation.
    double sum = 0.0;
    bool isPositive = true;
    for (int i = 1; i <= series; ++i) {
        if (isPositive)
            sum += std::pow((m - 1.0), i) / i;
        else
            sum -= std::pow((m - 1.0), i) / i;
        isPositive = !isPositive;
    }
    return sum + LN3 - LN2 + (j - 1023.0) * LN2;
}

/// @param series The series larger, result more precise
/// @note The m in range[1,2), √2m/2 in range[√2/2, √2). Let x moer close 1 for get higher precision.
double lnVersion2(double x, int series = 10) {
    if (x <= 0.0)
        throw std::domain_error("The x must is positive.");
    if (x == 1.0)
        return 0.0;
    unsigned long long _x = *reinterpret_cast<unsigned long long*>(&x);
    int j = _x >> 52;
    _x = 0x3ff0ULL << 48 ^ (_x << 12 >> 12);
    double m = *reinterpret_cast<double*>(&_x);
    m = std::pow(2.0, 0.5) * m / 2.0;
    double sum = 0.0;
    bool isPositive = true;
    for (int i = 1; i <= series; ++i) {
        if (isPositive)
            sum -= std::pow((m - 1.0), i) / i;
        else
            sum += std::pow((m - 1.0), i) / i;
        isPositive = !isPositive;
    }
    return sum + LN2 / 2.0 + (j - 1023.0) * LN2;
}

double log(double base, double x) {
    if (base == 1.0)
        throw std::domain_error("The base num can't is 1.");
    return lnVersion1(x) / lnVersion1(base);
}

}

namespace std_log
{

double ln(double x) {
    return std::log(x);
}

double log(double base, double x) {
    return ln(x) / ln(base);
}

}

namespace benchmark_test
{

void logPrecisionCompare(int testCount = 64, double start = 1.0, double step = 0.5) {
    std::cout.precision(12);

    for (int i = 0; i < testCount; ++i, start += step) {
        double a1 = my_log::lnVersion1(start);
        double a2 = std_log::ln(start);
        
        std::cout << "X: " << start << ";   \t";
        std::cout << "My Log: " << a1 << ";   \t";
        std::cout << "STD Log: " << a2 << ";   \t";
        std::cout << "Deviation: " << a2 - a1 << "\n";
    }
}

void logPerformanceCompare(int testCount = 1e6, double start = 1.0, double step = 0.5) {
    time_t startTime = clock();
    double x = start;
    for (int i = 0; i < testCount; ++i, x += step)
        std_log::ln(x);
    std::cout << "[STD][" << testCount << "]: ";
    std::cout << clock() - startTime << "msec" << '\n';
    
    startTime = clock();
    x = start;
    for (int i = 0; i < testCount; ++i, x += step)
        my_log::lnVersion1(x);
    std::cout << "[MY][" << testCount << "]: ";
    std::cout << clock() - startTime << "msec" << '\n';
}

}

int main()
{
    benchmark_test::logPrecisionCompare();
    benchmark_test::logPerformanceCompare(1e6);
    return 0;
}