<center><h1><b>在计算机程序中通过泰勒级数求解对数函数</b></h1></center>





## 一、自然对数函数的泰勒展开

$$
泰勒级数：\sum_{n=0}^{\infty}\frac{f^{(n)}(a)}{n!}(x-a)^n
\\\\newline
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

对于相同的自变量$$x$$，影响泰勒展开后结果精度的因素有二：

1. 泰勒展开的项数。展开项数越多，结果精度越高。
2. 自变量$$x$$与$$a$$（此处$$a=1$$）的差的绝对值大小（记作$$D$$）。$$D$$越小，即$$x$$与$$a$$值越接近，结果精度越高。

$$
对数函数的性质
\\
ln(a\cdot b)=ln(a)+ln(b)
\\
ln(a^b)=b\cdot ln(a)
$$

对于对数函数定义域$$(0,+\infty)$$内的任意$$x$$而言，使用上述性质将自变量标准化至1的邻域可以获得更精确的结果。
例如：在相同展开项数的情况下，直接对$$ln(2)$$进行泰勒展开求解的精准度不如$$4ln(2^{\frac{1}{4}})$$。

## 三、自变量为浮点数

对于浮点数而言，使用s代表其符号位，使用j代表其阶码，使用m代表其尾数。
由其存储方式可得单精度浮点数（下文默认使用单精度浮点数）
$$
F32=(-1)^s\times m\times\ 2^{(j-127)}
$$
或双精度浮点数
$$
F64=(-1)^s\times m\times\ 2^{(j-1023)}
$$
由于对数函数定义域为$$(0,+\infty)$$，故s始终为0，即
$$
F=m\times\ 2^{(j-127)}
$$
则浮点数作为自变量求解对数函数如下
$$
ln(F)=ln(m\times\ 2^{(j-127)})\\=ln(m)+ln(2^{(j-127)})\\=ln(m)+(j-127)\times ln(2)
$$
现需要将$$m\in[1,2)$$进行标准化，使其更加接近1，以增加计算精度。

<center><b>【法一】</b></center>

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
\ast 偏离量：0.6667
$$

<center><b>【法二】</b></center>

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
\ast 偏离量：0.7071
$$

## 四、ln(2)与ln(3)的求解

上文中用到了常数$$ln(2)$$与$$ln(3)$$，二者皆可通过将其标准化至1的领域进行泰勒展开求解。

在下面的代码实现中，函数*ln2ByTaylorSeries()*与*ln3ByTaylorSeries()*，通过以下方式将2与3标准化至1的领域。
$$
\S ln(2)
\\\\
\because ln(2^\frac{1}{16})=\frac{1}{16}ln(2)
\\
\therefore ln(2)=16ln(2^\frac{1}{16})
\\
\ast \quad 2-1=1, \ \sqrt[16]{2}-1\approx 0.0442737 \quad \ast

\\\\\\\\
\S ln(3)
\\\\
\because ln(3^\frac{1}{16})=\frac{1}{16}ln(3)
\\
\therefore ln(3)=16ln(3^\frac{1}{16})
\\
\ast \quad 3-1=2, \ \sqrt[16]{3}-1\approx 0.0710754 \quad \ast
$$

## 五、代码实现

```cpp
#include <cmath>
#include <exception>

namespace My
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
```

## 六、测试对比代码

```cpp
#include <iostream>

inline void logPrecisionCompare(int testCount = 20) {
    std::cout.precision(12);

    double x = 0.0;
    for (int i = 0; i < testCount; ++i) {
        x += 0.2;
        double a1 = MyLog::ln(x);
        double a2 = StdLog::ln(x);

        std::cout << "X: " << x << ";   \t";
        std::cout << "My ln: " << a1 << ";   \t";
        std::cout << "STD ln: " << a2 << ";   \t";
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
```

## 七、比对结果

![](./TestResults/2024-06-28_144900.png)

![](./TestResults/2024-06-28_144437.png)

1、速度上在十万的量级上有15ms左右的差距，标准库可能用到了常量表、优化算法、其他速度更快的计算方法。

2、在精准度上，MyLog的精准性是至少小数点后7位。
