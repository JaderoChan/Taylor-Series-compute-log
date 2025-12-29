#include <cmath>
#include <iostream>
#include <exception>

namespace my_log
{

template <typename T>
constexpr double powInt(T base, int exp) {
    if (exp == 0)   return 1.0;
    if (exp < 0)    return 1.0 / powInt(base, -exp);

    double result = 1.0;
    while (exp) {
        if (exp & 1)
            result *= base;
        base *= base;
        exp >>= 1;
    }

    return result;
}

constexpr double LN2 = 0.693147180559945309417232121458;
constexpr double LN3 = 1.098612288668109691395245236922;
constexpr double SQRT2 = 1.414213562373095048801688724209;

/// @param series The series larger, result more precise.
/// @note The m in range[1,2), 2m/3 in range[2/3, 4/3). Let "x" moer close 1 for higher precision.
double lnVersion1(double x, int series = 10) {
    if (x <= 0.0)   // The log() domain of definition is (0, +infty).
        throw std::domain_error("The x must is positive.");
    if (x == 1.0)   // The ln(1) = 0.
        return 0.0;
    // Cast the double to unsigned long long for bit operation.
    unsigned long long v;
    static_assert(sizeof(v) == sizeof(x));
    memcpy(&v, &x, sizeof(v));
    // Get the float point num's exp.
    int j = v >> 52;
    // Get the float point num's mantissa. The (0x3ff0ULL << 48) is a double it = 1.0.
    v = 0x3ff0ULL << 48 ^ (v << 12 >> 12);
    // Cast the bit operated unsigned long long to double.
    double m = *reinterpret_cast<double*>(&v);
    // Let "x" more close 1 for higher precision.
    // m = 2.0 * m / 3.0;
    m *= 0.666666666666666666666666666666;
    // The Taylor series expanded summation.
    double result = 0.0;
    int sign = 1;
    for (int i = 1; i <= series; ++i) {
        result += sign * powInt((m - 1.0), i) / i;
        sign *= -1;
    }
    return result + LN3 - LN2 + (double) (j - 1023) * LN2;
}

/// @param series The series larger, result more precise
/// @note The m in range[1,2), √2m/2 in range[√2/2, √2). Let x moer close 1 for get higher precision.
double lnVersion2(double x, int series = 10) {
    if (x <= 0.0)
        throw std::domain_error("The x must is positive.");
    if (x == 1.0)
        return 0.0;
    unsigned long long v;
    static_assert(sizeof(v) == sizeof(x));
    memcpy(&v, &x, sizeof(v));
    int j = v >> 52;
    v = 0x3ff0ULL << 48 ^ (v << 12 >> 12);
    double m = *reinterpret_cast<double*>(&v);
    m = SQRT2 * m / 2.0;
    double result = 0.0;
    int sign = 1;
    for (int i = 1; i <= series; ++i) {
        result += sign * powInt((m - 1.0), i) / i;
        sign *= -1;
    }
    return result + LN2 / 2.0 + (double) (j - 1023) * LN2;
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

namespace benchmark
{

void lnPrecisionCompare(int testCount = 64, double start = 1.0, double step = 0.5) {
    std::cout.precision(12);

    for (int i = 0; i < testCount; ++i, start += step) {
        double a1 = my_log::lnVersion1(start, 1000);
        double a2 = std_log::ln(start);

        std::cout << "X: " << start << ";   \t";
        std::cout << "My Log: " << a1 << ";   \t";
        std::cout << "STD Log: " << a2 << ";   \t";
        std::cout << "Deviation: " << a2 - a1 << "\n";
    }
}

void lnPerformanceCompare(int testCount = 1e6, double start = 1.0, double step = 0.5) {
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
    benchmark::lnPrecisionCompare();
    benchmark::lnPerformanceCompare(1e6);
    return 0;
}
