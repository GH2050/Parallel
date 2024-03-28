#include <iostream>
#include <Windows.h>// 可更换成Linux的计时方法gettimeofday()
#include<iomanip>

const long long Nsize = 1048576;// 1024*1024; // 2^13 2^16 2^18 double for 64KB(128KB) 2*256KB 3MB // float则数量乘二
float vec[Nsize]{};
float sum = 0;
float sum1 = 0, sum2 = 0;

void Initial(long long Nsize) {
    for (long long i = 0; i < Nsize; i++) {
        switch (i % 3)
        {
        case 0:
            vec[i] = static_cast<float>(1) / 3;
            break;
        case 1:
            vec[i] = static_cast<float>(1) / 7;
            break;
        case 2:
            vec[i] = static_cast<float>(1) / 11;
            break;
        }
    }

}

// 平凡算法：累加求和
void Simple(long long Nsize) {
    for (long long i = 0; i < Nsize; i++) {
        sum += vec[i];
    }
}
// 平凡算法：累加求和
void Simple2(long long Nsize) {
    for (long long i = 0; i < Nsize; i += 2) {
        sum += vec[i];
        sum += vec[i + 1];
    }
}

// 多路链式：流水线
void Pipline(long long Nsize){
    for (long long i = 0; i < Nsize;i += 2) {
        sum1 += vec[i];
        sum2 += vec[i + 1];
    }
    sum = sum1 + sum2;
}

// 循环展开
template<int N>
struct LoopUnroller {
    static void unroll() {
        sum += vec[Nsize - N];
        LoopUnroller<N - 1>::unroll();
    }
};
template<>
struct LoopUnroller<1> {
    static void unroll() {
        ;
    }
};

// 递归循环
void Recursion(long long Nsize) {
    if (Nsize == 1) {
        sum = vec[0];
        return;
    }
    else
    {
        for (long long i = 0; i < Nsize; i++)
            vec[i] += vec[Nsize - i - 1];
        Nsize /= 2;
        Recursion(Nsize);
    }
}

// 二重循环(Recursion)
void NoRecursion(long long Nsize) {
    for (long long i = Nsize; i > 1; i /= 2)
        for (long long j = 0; j < i / 2; j++)
            vec[j] += vec[i - j - 1];
    sum = vec[0];
}

// 二重循环2
void NoRecursion2(long long Nsize) {
    for (long long i = Nsize; i > 1; i /= 2)
        for (long long j = 0; j < i / 2; j++)
            vec[j] = vec[j * 2] + vec[j * 2 + 1];
    sum = vec[0];
}

// 二重循环3(Simple)
void NoRecursion3(long long Nsize) {
    for (long long i = Nsize; i > 1; i /= 2)
        for (long long j = 0; j < i / 2; j++)
            sum += vec[Nsize - i + j];
    sum += vec[Nsize - 1];
}

// Test01 测试数据规模对性能的影响 // 除了模板展开，因为模板展开需要提前确定大小
void Size_Test() {
    // 测试数据规模变化对程序运行的影响
    // 多种方法的平均运行时间

    long long min_test_size = 1024;// 1024
    long long max_test_size = 1024;
    max_test_size *= 1024;
    max_test_size *= 1;
    int test_step = 2;// 每次规模翻倍
    LARGE_INTEGER frequency, t1, t2;
    long long repeat = 0;// 用来对运行时间较短的程序进行循环运行
    double elapsedTime = 0;

    // 获取时钟频率
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // 进行测试并输出结果
    for (long long n = min_test_size; n <= max_test_size; n *= test_step)
    {
        // 进行测试 平凡算法
        Initial(max_test_size);
        elapsedTime = 0;
        for (repeat = 1; ; repeat++)// Simple
        {
            //Initial(max_test_size);
            sum = 0; sum1 = 0; sum2 = 0;
            QueryPerformanceCounter(&t1);
            Simple(n);
            QueryPerformanceCounter(&t2);
            elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
            if (repeat >= 1000 || elapsedTime > 100)
            {
                break;
            }
        }
        elapsedTime /= repeat; // 单位为毫秒
        std::cout << std::fixed << std::setprecision(20) << std::setw(15) << "Simple," << std::setw(8) << n << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
        std::cout << sum << std::endl;// << std::fixed << std::setprecision(10) << 永久效果

        // 进行测试 平凡算法2
        Initial(max_test_size);
        elapsedTime = 0;
        for (repeat = 1; ; repeat++)// Simple2
        {
            //Initial(max_test_size);
            sum = 0; sum1 = 0; sum2 = 0;
            QueryPerformanceCounter(&t1);
            Simple2(n);
            QueryPerformanceCounter(&t2);
            elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
            if (repeat >= 1000 || elapsedTime > 100)
            {
                break;
            }
        }
        elapsedTime /= repeat; // 单位为毫秒
        std::cout << std::setw(15) << "Simple2," << std::setw(8) << n << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
        std::cout << sum << std::endl;

        // 进行测试 多路链式
        elapsedTime = 0;
        for (repeat = 1; ; repeat++)// Pipline
        {
            //Initial(max_test_size);
            sum = 0; sum1 = 0; sum2 = 0;
            QueryPerformanceCounter(&t1);
            Pipline(n);
            QueryPerformanceCounter(&t2);
            elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
            if (repeat >= 1000 || elapsedTime > 100)
            {
                break;
            }
        }
        elapsedTime /= repeat; // 单位为毫秒
        std::cout << std::setw(15) << "Pipline," << std::setw(8) << n << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
        std::cout << sum  << std::endl;

        // 进行测试 递归循环
        elapsedTime = 0;
        for (repeat = 1; ; repeat++)// Recursion
        {
            Initial(max_test_size);
            sum = 0; sum1 = 0; sum2 = 0;
            QueryPerformanceCounter(&t1);
            Recursion(n);
            QueryPerformanceCounter(&t2);
            elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
            if (repeat >= 1000 || elapsedTime > 100)
            {
                break;
            }
        }
        elapsedTime /= repeat; // 单位为毫秒
        std::cout << std::setw(15) << "Recursion," << std::setw(8) << n << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
        std::cout << sum  << std::endl;

        // 进行测试 二重循环
        elapsedTime = 0;
        for (repeat = 1; ; repeat++)// NoRecursion
        {
            Initial(max_test_size);
            sum = 0; sum1 = 0; sum2 = 0;
            QueryPerformanceCounter(&t1);
            NoRecursion(n);
            QueryPerformanceCounter(&t2);
            elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
            if (repeat >= 1000 || elapsedTime > 100)
            {
                break;
            }
        }
        elapsedTime /= repeat; // 单位为毫秒
        std::cout << std::setw(15) << "NoRecursion," << std::setw(8) << n << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
        std::cout << sum  << std::endl;

        // 进行测试 二重循环2
        elapsedTime = 0;
        for (repeat = 1; ; repeat++)// NoRecursion2
        {
            Initial(max_test_size);
            sum = 0; sum1 = 0; sum2 = 0;
            QueryPerformanceCounter(&t1);
            NoRecursion2(n);
            QueryPerformanceCounter(&t2);
            elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
            if (repeat >= 1000 || elapsedTime > 100)
            {
                break;
            }
        }
        elapsedTime /= repeat; // 单位为毫秒
        std::cout << std::setw(15) << "NoRecursion2," << std::setw(8) << n << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
        std::cout << sum << std::endl;

        // 进行测试 二重循环3
        elapsedTime = 0;
        for (repeat = 1; ; repeat++)// NoRecursion3
        {
            Initial(max_test_size);
            sum = 0; sum1 = 0; sum2 = 0;
            QueryPerformanceCounter(&t1);
            NoRecursion3(n);
            QueryPerformanceCounter(&t2);
            elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
            if (repeat >= 1000 || elapsedTime > 100)
            {
                break;
            }
        }
        elapsedTime /= repeat; // 单位为毫秒
        std::cout << std::setw(15) << "NoRecursion3," << std::setw(8) << n << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
        std::cout << sum << std::endl;
    }

}

// Test02 测试模板展开，因为模板展开需要提前确定大小 // 编译需要选项-ftemplate-depth=
void Size_Test2() {
    // 测试数据规模变化对程序运行的影响

    LARGE_INTEGER frequency, t1, t2;
    long long repeat = 0;// 用来对运行时间较短的程序进行循环运行
    double elapsedTime = 0;

    // 获取时钟频率
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    //LoopUnroller<1024>::unroll();
    //LoopUnroller<2048>::unroll();
    //LoopUnroller<4096>::unroll();
    //LoopUnroller<8192>::unroll();
    //LoopUnroller<16384>::unroll();
    //LoopUnroller<32768>::unroll();
    //LoopUnroller<65536>::unroll();
    //LoopUnroller<131072>::unroll();
    //LoopUnroller<262144>::unroll();
    //LoopUnroller<524288>::unroll();
    //LoopUnroller<1048576>::unroll();

    // 进行测试
    Initial(1024);
    elapsedTime = 0;
    for (repeat = 1; ; repeat++)
    {
        sum = 0; sum1 = 0; sum2 = 0;
        QueryPerformanceCounter(&t1);
        LoopUnroller<1024>::unroll();
        QueryPerformanceCounter(&t2);
        elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
        if (repeat >= 1000 || elapsedTime > 100)
        {
            break;
        }
    }
    elapsedTime /= repeat; // 单位为毫秒
    std::cout << std::setw(15) << "LoopUnroller," << std::setw(8) << 1024 << ",循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
    std::cout << sum << std::endl;

}

// Test03 对1024进行测试，尝试发现L3_Miss成为主流 // 专注单一程序单一规模
void Perf_Test() {

    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// 用来对运行时间较短的程序进行循环运行
    double elapsedTime = 0;

    // 获取时钟频率
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // 进行测试 模板展开
//    Initial(1024);
//    elapsedTime = 0;
//    for (repeat = 1; ; repeat++)
//    {
//        sum = 0; sum1 = 0; sum2 = 0;
//        QueryPerformanceCounter(&t1);
//        LoopUnroller<1024>::unroll();
//        QueryPerformanceCounter(&t2);
//        elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
//        if (repeat >= 1000 || elapsedTime > 100)
//        {
//            break;
//        }
//    }
//    elapsedTime /= repeat; // 单位为毫秒
//    std::cout << std::setw(15) << "LoopUnroller," << "循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
//    std::cout << sum << std::endl;

//    // 进行测试 平凡算法
//    Initial(1024);
//    elapsedTime = 0;
//    for (repeat = 1; ; repeat++)// Simple
//    {
//        //Initial(1024);
//        sum = 0; sum1 = 0; sum2 = 0;
//        QueryPerformanceCounter(&t1);
//        Simple(1024);
//        QueryPerformanceCounter(&t2);
//        elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
//        if (repeat >= 1000 || elapsedTime > 100)
//        {
//            break;
//        }
//    }
//    elapsedTime /= repeat; // 单位为毫秒
//    std::cout << std::setw(15) << "Simple," << "循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
//    std::cout << sum << std::endl;
//
    // 进行测试 递归循环
//    elapsedTime = 0;
//    for (repeat = 1; ; repeat++)// Recursion
//    {
//        Initial(1024);
//        sum = 0; sum1 = 0; sum2 = 0;
//        QueryPerformanceCounter(&t1);
//        Recursion(1024);
//        QueryPerformanceCounter(&t2);
//        elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
//        if (repeat >= 1000 || elapsedTime > 100)
//        {
//            break;
//        }
//    }
//    elapsedTime /= repeat; // 单位为毫秒
//    std::cout << std::setw(15) << "Recursion," << "循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
//    std::cout << sum << std::endl;

    // 进行测试 多路链式
    elapsedTime = 0;
    for (repeat = 1; ; repeat++)// Recursion
    {
        Initial(1024);
        sum = 0; sum1 = 0; sum2 = 0;
        QueryPerformanceCounter(&t1);
        Pipline(1024);
        QueryPerformanceCounter(&t2);
        elapsedTime += (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
        if (repeat >= 1000 || elapsedTime > 100)
        {
            break;
        }
    }
    elapsedTime /= repeat; // 单位为毫秒
    std::cout << std::setw(15) << "Pipline," << "循环次数," << std::setw(6) << repeat << ",平均用时," << std::setw(30) << elapsedTime << ",ms,";
    std::cout << sum << std::endl;

}

int main()
{
    Initial(Nsize);
    //Size_Test();
    //Size_Test2();
    Perf_Test();
    return 0;
}
