#include <iostream>
#include <Windows.h>
#include<iomanip>

const int Nsize = 2000;// 10000;//120 500 1400 for 64KB(128KB) 2*256KB 3MB
int vec[Nsize]{};
int mat[Nsize][Nsize]{};
int sum[Nsize]{};

void Initial(int Nsize) {
    for (int i = 0; i < Nsize; i++)
    {
        for (int j = i; j < Nsize; j++)
        {
            mat[i][j] = 1;
        }
    }
    for (int i = 0; i < Nsize; i++)
    {
        vec[i] = 1;
    }
    for (int j = 0; j < Nsize; j++)
    {
        sum[j] = 0;
    }
    // check 矩阵是否正确初始化
    //std::cout << "初始矩阵:" << std::endl;
    //for (int i = 0; i < Nsize; i++)
    //{
    //    for (int j = 0; j < Nsize; j++) {
    //        std::cout << mat[i][j];
    //        if (j < Nsize - 1)
    //            std::cout << ' ';
    //        else
    //            std::cout << std::endl;
    //    }
    //}
}

// 清零计算结果
void Zero_Sum(int Nsize) {
    for (int j = 0; j < Nsize; j++)
    {
        sum[j] = 0;
    }
}

// 输出计算结果
void Check_Sum(int Nsize, int repeat) {
    for (int i = 0; i < Nsize; i++)
    {
        std::cout << sum[i] / repeat;
        if (i < Nsize - 1)
            std::cout << ' ';
        else
            std::cout << std::endl;
    }
}

// 平凡算法：按列求和
void Col_Calc(int Nsize) {
    for (int j = 0; j < Nsize; j++)
    {
        for (int i = 0; i < Nsize; i++)
        {
            sum[j] += vec[i] * mat[i][j];
        }
    }
}

// Cache算法：按行求和
void Row_Calc(int Nsize) {
    for (int i = 0; i < Nsize; i++)
    {
        for (int j = 0; j < Nsize; j++)
        {
            sum[j] += vec[i] * mat[i][j];
        }
    }
}



// Test01 测试数据规模对性能的影响
void Size_Test() {
    // 测试数据规模变化对程序运行的影响
    // 两种方法的平均运行时间

    int min_test_size = 10;// 100 -> 10000 // 10 -> 2000
    int max_test_size = 2000;
    int test_step = 10;// 50;// 不改变步长
    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// 用来对运行时间较短的程序进行循环运行
    double elapsedTime = 0;

    // 获取时钟频率
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // 进行测试并输出结果
    for (int n = min_test_size; n <= max_test_size; n+=test_step)
    {
        // 进行测试 列计算
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Col_Calc
        {
            Col_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
        std::cout << "列计算 规模: " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";

        // 进行测试 行计算
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Row_Calc
        {
            Row_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
        std::cout << "行计算 规模: " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";
    }

}

// Test02 测试CacheLine大小，通过手动调整Row_Calc的大小
void CacheLine_Size_Test() {

    int min_test_size = 2;// 2 -> 200
    int max_test_size = 200;
    int test_step = 2;// 不调整步长
    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// 用来对运行时间较短的程序进行循环运行
    double elapsedTime = 0;

    // 获取时钟频率
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // 进行测试并输出结果
    for (int n = min_test_size; n <= max_test_size; n += test_step)
    {
        // 进行测试 列计算
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Col_Calc
        {
            Col_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
        std::cout << "列计算 规模: " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";

        // 进行测试 行计算
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Row_Calc
        {
            Row_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
        std::cout << "行计算 规模: " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";

        // 不调整步长
    }

}

// Test03 利用Vtune进行profile前 200-800 与 800-1400 在500与1100附近进行测试 重复CacheLine测试方法 尝试证明明显性能下降与L2与L3 Cache命中率相关
void Pre_Vtune_Test() {

    int min_test_size = 800;// 800 -> 1400 // 200 -> 800
    int max_test_size = 1400;
    int test_step = 8;// 不调整步长
    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// 用来对运行时间较短的程序进行循环运行
    double elapsedTime = 0;

    // 获取时钟频率
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // 进行测试并输出结果
    for (int n = min_test_size; n <= max_test_size; n += test_step)
    {
        // 进行测试 列计算
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Col_Calc
        {
            Col_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
        std::cout << "Col_Calc " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";

        // 进行测试 行计算
        Zero_Sum(max_test_size);
        QueryPerformanceCounter(&t1);
        for (repeat = 1; ; repeat++)// Row_Calc
        {
            Row_Calc(n);
            QueryPerformanceCounter(&t2);
            if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
            {
                break;
            }
        }
        elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
        std::cout << "Row_Calc " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";

        // 不调整步长
    }

}

// Test04 对900 1100 1300进行测试，尝试发现L3_Miss成为主流 // 专注单一程序单一规模
void Vtune_Test() {

    int min_test_size = 900;// 900 -> 1300
    int max_test_size = 1300;// 选取合适的三个点位
    int test_step = 200;// 不调整步长
    int n = min_test_size + test_step * 1;// 修改规模 注释掉其中一个程序

    LARGE_INTEGER frequency, t1, t2;
    int repeat = 0;// 用来对运行时间较短的程序进行循环运行
    double elapsedTime = 0;

    // 获取时钟频率
    if (!QueryPerformanceFrequency(&frequency)) {
        std::cout << "QueryPerformanceFrequency failed!\n";
        return;
    }

    // 进行测试 列计算
    Zero_Sum(max_test_size);
    QueryPerformanceCounter(&t1);
    for (repeat = 1; ; repeat++)// Col_Calc
    {
        Col_Calc(n);
        QueryPerformanceCounter(&t2);
        if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
        {
            break;
        }
    }
    elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
    std::cout << "列计算 规模: " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";

    // 进行测试 行计算
//    Zero_Sum(max_test_size);
//    QueryPerformanceCounter(&t1);
//    for (repeat = 1; ; repeat++)// Row_Calc
//    {
//        Row_Calc(n);
//        QueryPerformanceCounter(&t2);
//        if ((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart > 100)// 循环100ms
//        {
//            break;
//        }
//    }
//    elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart / repeat; // 单位为毫秒
//    std::cout << "行计算 规模: " << std::setw(4) << n << " 循环次数: " << std::setw(6) << repeat << " 平均用时: " << elapsedTime << " ms\n";

}

int main()
{
    // Initial Nsize = 2000;
    Initial(Nsize);

    // 测试数据 规模对性能的影响 分析Cache大小 引出Cache Line问题
    // Size_Test();

    // 探索Cache Line大小(N = 2 -> 200)
    // CacheLine_Size_Test();

    // 利用Vtune进行profile前 在1100附近(N = 800 -> 1400)进行测试 重复CacheLine测试方法 尝试证明明显性能下降与L3 Cache命中率相关
    // Pre_Vtune_Test();

    // 结合Vtune分析两者在Cache命中率上的差异 证明命中率及L3大小极大地影响了运行速度，即主存取影响了速度
    Vtune_Test();


    return 0;
}
