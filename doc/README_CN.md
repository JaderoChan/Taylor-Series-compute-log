# 在计算机程序中通过泰勒级数求解对数函数

[**简体中文** | [**English**](./README_EN.md)]

## 介绍

不借助任何第三方库与标准数学库，通过泰勒级数对对数函数进行求解。

## 一、自然对数函数的泰勒展开

$$
泰勒级数：\sum_{n=0}^{\infty}\frac{f^{(n)}(a)}{n!}(x-a)^n
\\\\
对f(x)=ln(x)\quad(a=1)进行泰勒展开得
\\\\
\Rightarrow  \frac{f^{(0)}(1)}{0!}(x-1)^0+\frac{f^{(1)}(1)}{1!}(x-1)^1+\frac{f^{(2)}(1)}{2!}(x-1)^2+\cdots+\frac{f^{(n)}(1)}{n!}(x-1)^n
\\\\
\Rightarrow  0+(x-1)+(-\frac{1}{2})(x-1)^2+(\frac{2}{2\times3})(x-1)^3+(-\frac{6}{2\times3\times4})(x-1)^4+\cdots
\\\\
\Rightarrow  \frac{(x-1)^1}{1}-\frac{(x-1)^2}{2}+\frac{(x-1)^3}{3}-\frac{(x-1)^4}{4}+\cdots
\\\\
\Rightarrow  \sum_{n=1}^{\infty}(-1)^{n-1}\frac{(x-1)^n}{n}
$$

## 二、有关自然对数函数泰勒展开求解的精度探讨

对于相同的自变量$x$，影响泰勒展开后结果精度的因素有二：

1. 泰勒展开的项数。展开项数越多，结果精度越高。
2. 自变量$x$与$a$（此处$a=1$）的差的绝对值大小（记作$D$）。$D$越小，即$x$与$a$值越接近，结果精度越高。

$$
对数函数的性质
\\
ln(a\cdot b)=ln(a)+ln(b)
\\
ln(a^b)=b\cdot ln(a)
$$

对于对数函数定义域$(0,+\infty)$内的任意$x$而言，使用上述性质将自变量标准化至1的邻域可以获得更精确的结果。
例如：在相同展开项数的情况下，直接对$ln(2)$进行泰勒展开求解的精准度不如$4ln(2^{\frac{1}{4}})$。

## 三、自变量为浮点数

对于浮点数而言，使用s代表其符号位，使用j代表其阶码，使用m代表其尾数。
由其存储方式可得单精度浮点数（下文默认使用单精度浮点数）
$$
F32=(-1)^s\times m\times 2^{(j-127)}
$$
或双精度浮点数
$$
F64=(-1)^s\times m\times 2^{(j-1023)}
$$
由于对数函数定义域为$(0,+\infty)$，故s始终为0，即
$$
F=m\times 2^{(j-127)}
$$
则浮点数作为自变量求解对数函数如下
$$
ln(F)=ln(m\times 2^{(j-127)})
\\
=ln(m)+ln(2^{(j-127)})
\\
=ln(m)+(j-127)\times ln(2)
$$
现需要将$m\in[1,2)$进行标准化，使其更加接近1，以增加计算精度。

<p><b><center>【法一】</center></b></p>

$$
ln(\frac{3}{2}\cdot\frac{2}{3}\cdot m)=ln(3)-ln(2)+ln(\frac{2}{3}m)
\\\\
\because m\in[1,2)
\\\\
\therefore \frac{2}{3}m\in\left([\frac{2}{3},\frac{4}{3})\approx[0.666,1.333)\right)
\\\\
ln(F)=ln(3)-ln(2)+(j-127)\times ln(2)+ln(\frac{2}{3}m)
\\
其中ln(2)，ln(3)皆为常数，(j-127)为浮点数阶码偏移值。
\\
故只需对ln(\frac{2}{3}m)进行泰勒展开求和即可
\\\\
\ast 偏离量：0.6667 \ast
$$

<p><b><center>【法二】</center></b></p>

$$
ln(\sqrt{2}\cdot\frac{\sqrt{2}}{2}\cdot m)=\frac{1}{2}ln(2)+ln(\frac{\sqrt{2}}{2}m)
\\\\
\because m\in[1,2)
\\\\
\therefore \frac{\sqrt{2}}{2}m\in\left([\frac{\sqrt{2}}{2},\sqrt{2})\approx[0.707,1.414)\right)
\\\\
ln(F)=\frac{1}{2}ln(2)+(j-127)\times ln(2)+ln(\frac{\sqrt{2}}{2}m)
\\
其中ln(2)为常数，(j-127)为浮点数阶码偏移值。
\\
故只需对ln(\frac{\sqrt{2}}{2}m)进行泰勒展开求和即可
\\\\
\ast 偏离量：0.7071 \ast
$$

## 四、ln(2)与ln(3)的求解

上文中用到了常数$ln(2)$与$ln(3)$，二者皆可通过将其标准化至1的领域进行泰勒展开求解。

在下面的代码实现中，函数*ln2ByTaylorSeries()*与*ln3ByTaylorSeries()*，通过以下方式将2与3标准化至1的邻域。
$$
--ln(2)--
\\\\
\because ln(2^\frac{1}{16})=\frac{1}{16}ln(2)
\\
\therefore ln(2)=16ln(2^\frac{1}{16})
\\
\ast \quad 2-1=1, \quad \sqrt[16]{2}-1\approx 0.0442737 \quad \ast
\\\\
--ln(3)--
\\\\
\because ln(3^\frac{1}{16})=\frac{1}{16}ln(3)
\\
\therefore ln(3)=16ln(3^\frac{1}{16})
\\
\ast \quad 3-1=2, \quad \sqrt[16]{3}-1\approx 0.0710754 \quad \ast
$$

## 五、代码实现

```c
#include <assert.h> // assert()
#include <string.h> // memcpy()

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
/// @note The m in range[1,2), 2m/3 in range[2/3, 4/3). Let "x" (see also above document) moer close 1.0 for higher precision.
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
/// @note The m in range[1,2), √2m/2 in range[√2/2, √2). Let "x" (see also above document) moer close 1.0 for higher precision.
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
```

## 六、测试对比代码

```c
#include <math.h>   // log() for contrast difference.
#include <time.h>   // time_t, clock()
#include <stdio.h>  // printf()

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
```

## 七、比对结果

![benchmark](../benchmark_image/benchmark.png)
