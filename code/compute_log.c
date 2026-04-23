#include <assert.h> // assert()
#include <math.h>   // log() for contrast difference.
#include <time.h>   // time_t, clock()
#include <string.h> // memcpy()
#include <stdio.h>  // printf()

#define ABS(x) ((x) < 0.0 ? -(x) : (x))

// some precomputed constants.
#define LN2     0.693147180559945309417232121458
#define LN3     1.098612288668109691395245236922
#define SQRT2   1.414213562373095048801688724209

/// @brief Compute integer exponential powers using fast exponentiation.
double my_pow(double x, int exp)
{
    if (exp == 0)
        return 1.0;
    if (exp < 0)
        return 1.0 / my_pow(x, -exp);

    double result = 1.0;
    while (exp)
    {
        if (exp & 1)
            result *= x;
        x *= x;
        exp >>= 1;
    }

    return result;
}

/// @param series The series larger, result more precise.
/// @note The m in range[1,2), 2m/3 in range[2/3, 4/3). See also document.
double my_ln_ver1(double x, int series)
{
    // The ln() domain of definition is (0, +infty).
    assert(x > 0.0);

    // ln(1.0) = 0.0.
    if (x == 1.0)
        return 0.0;

    // Bit-level reinterpretation of double type into unsigned long long type (similar to reinterpret_cast in C++).
    unsigned long long v;
    memcpy(&v, &x, sizeof(v));

    // Get the float point num's exp.
    int j = v >> 52;
    // Get the float point num's mantissa. The (0x3ff0ULL << 48) is a double, it is 1.0.
    v = 0x3ff0ULL << 48 ^ (v << 12 >> 12);

    // Bit-level reinterpretation of unsigned long long type into double type.
    memcpy(&x, &v, sizeof(x));

    // Let "x" more close 1.0 for higher precision.
    // x = 2.0 * x / 3.0;
    x *= 0.666666666666666666666666666666;

    // The Taylor series expanded summation.
    double result = 0.0;
    int sign = 1;
    for (int i = 1; i <= series; ++i)
    {
        result += sign * my_pow((x - 1.0), i) / i;
        sign *= -1;
    }

    return result + LN3 - LN2 + (double) (j - 1023) * LN2;
}

/// @param series The series larger, result more precise.
/// @note The m in range[1,2), √2m/2 in range[√2/2, √2). See also document.
double my_ln_ver2(double x, int series)
{
    assert(x > 0.0);

    if (x == 1.0)
        return 0.0;

    unsigned long long v;
    memcpy(&v, &x, sizeof(v));

    int j = v >> 52;
    v = 0x3ff0ULL << 48 ^ (v << 12 >> 12);

    memcpy(&x, &v, sizeof(x));

    // Let "x" more close 1.0 for higher precision.
    // x = SQRT2 * x / 2.0;
    x *= 0.70710678118654752440084436210485;

    // The Taylor series expanded summation.
    double result = 0.0;
    int sign = 1;
    for (int i = 1; i <= series; ++i)
    {
        result += sign * my_pow((x - 1.0), i) / i;
        sign *= -1;
    }

    return result + LN2 / 2.0 + (double) (j - 1023) * LN2;
}

double std_ln(double x)
{
    return log(x);
}

void benchmark_precision_compare(
    double start, double end, double step,
    double (*my_ln)(double, int), int my_ln_series)
{
    for (; start <= end; start += step)
    {
        double my_ans = my_ln(start, my_ln_series);
        double std_ans = std_ln(start);

        printf("X: %lf;\tMy ln(): %lf;\tStd ln(): %lf;\tDeviation: %lf.\n",
            start, my_ans, std_ans, ABS(my_ans - std_ans));
    }
}

void benchmark_performance_compare(
    double start, double end, double step,
    double (*my_ln)(double, int), int my_ln_series)
{
    int count = (int) ((end - start) / step) + 2;

    // My ln()
    time_t start_time = clock();

    double x = start;
    for (; x <= end; x += step)
        my_ln(x, my_ln_series);

    time_t end_time = clock();
    printf("[My ln()][%d] %ld msec.\n", count, end_time - start_time);

    // Std ln()
    start_time = clock();

    x = start;
    for (; x <= end; x += step)
        std_ln(x);

    end_time = clock();
    printf("[Std ln()][%d] %ld msec.\n", count, end_time - start_time);
}

int main(int argc, char* argv[])
{
    // Ensure that the double type can be correctly converted to unsigned long long
    // to extract the exponent, mantissa, etc. of the double.
    assert(sizeof(unsigned long long) == sizeof(double));

    benchmark_precision_compare(1.0, 100.0, 0.01, my_ln_ver1, 100);
    benchmark_performance_compare(1.0, 10000.0, 0.5, my_ln_ver1, 100);

    return 0;
}
