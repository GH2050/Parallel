// 在Windows、x86下进行测试,不支持aligned_alloc但是支持_aligned_malloc
// 在另一文件中使用Linux、arm进行测试，可以更改数据规模并实现对齐的内存分配
#include <iostream>

// Windows
#include <windows.h> // 传统的min/max宏定义会影响程序中使用std::min/max，可以使用(std::min/max)阻止宏替换
#include <cstdlib> // _aligned
// Linux
//#include <time.h> // gettimeofday

// arm
//#include <arm_neon.h> // NEON

// x86
#include <nmmintrin.h> // SSE4.2
#include <immintrin.h> // AVX2


void Generate(int* Weight, int* Value, int N) {
    for (int i = 0; i < N; i++) {
        Weight[i] = rand() % 100;
        Value[i] = rand() % 100;
    }
}

void Zero(int* Array, int N) {
    for (int i = 0; i < N; i++) {
        Array[i] = 0;
    }
}

// 串行算法
void Knapsack_Basic(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        for(int j = 0; j <= Capacity; j++) {
            if (j < Weight[i]) {
                Ptr_sol_cur[j] = -Value[i];
            }
            else {
                Ptr_sol_cur[j] = Ptr_sol_pre[j-Weight[i]];// 也许可以并行移动
            }
        }

        for(int j = 0; j <= Capacity; j++) {// 可以进行并行优化
            Ptr_sol_cur[j] += Value[i];
            Ptr_sol_cur[j] = (std::max)(Ptr_sol_pre[j], Ptr_sol_cur[j]);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;

    }

    // 输出结果
    //std::cout << "Basic Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// 串行算法省略不需要的计算
void Knapsack_Basic_update(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {

        for(int j = Weight[i]; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-Weight[i]];
        }

        for(int j = 0; j < Weight[i]; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        for(int j = Weight[i]; j <= Capacity; j++) {
            Ptr_sol_cur[j] += Value[i];
            Ptr_sol_cur[j] = (std::max)(Ptr_sol_pre[j], Ptr_sol_cur[j]);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "Basic Update Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// 串行算法提前设定临时变量存储Value[i]和Weight[i]避免运行中频繁访问数组
void Knapsack_Basic_update_temp(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    int val, wei;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        val = Value[i];
        wei = Weight[i];
        for(int j = wei; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-wei];
        }

        for(int j = 0; j < wei; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        for(int j = wei; j <= Capacity; j++) {
            Ptr_sol_cur[j] += val;
            Ptr_sol_cur[j] = (std::max)(Ptr_sol_pre[j], Ptr_sol_cur[j]);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "Basic Update Temp Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// SSE并行算法
void Knapsack_SSE(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m128i sol_pre, sol_cur, val;


    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        for(int j = 0; j <= Capacity; j++) {
            if (j < Weight[i]) {
                Ptr_sol_cur[j] = -Value[i];
            }
            else {
                Ptr_sol_cur[j] = Ptr_sol_pre[j-Weight[i]];
            }
        }

        // 指定val为Value[i]
        val = _mm_set1_epi32(Value[i]);
        for(int j = 0; j <= Capacity; j+=4) {
            sol_cur = _mm_loadu_si128((__m128i*)& Ptr_sol_cur[j]);
            //val = _mm_loadu_si128((__m128i*)&Value[i]);
            sol_cur = _mm_add_epi32(sol_cur, val);

            sol_pre = _mm_loadu_si128((__m128i*)& Ptr_sol_pre[j]);
            sol_cur = _mm_max_epi32(sol_pre, sol_cur);
            _mm_storeu_si128((__m128i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "SSE Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// SSE并行算法省略不需要的计算
void Knapsack_SSE_update(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m128i sol_pre, sol_cur, val;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        for(int j = Weight[i]; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-Weight[i]];
        }

        for(int j = 0; j < Weight[i]; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        // 指定val为Value[i]
        val = _mm_set1_epi32(Value[i]);
        for(int j = Weight[i]; j <= Capacity; j+=4) {
            sol_cur = _mm_loadu_si128((__m128i*)& Ptr_sol_cur[j]);
            //val = _mm_loadu_si128((__m128i*)&Value[i]);
            sol_cur = _mm_add_epi32(sol_cur, val);

            sol_pre = _mm_loadu_si128((__m128i*)& Ptr_sol_pre[j]);
            sol_cur = _mm_max_epi32(sol_pre, sol_cur);
            _mm_storeu_si128((__m128i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "SSE Update Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// SSE并行算法提前设定临时变量存储Value[i]和Weight[i]避免运行中频繁访问数组
void Knapsack_SSE_update_temp(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m128i sol_pre, sol_cur, val;
    int val_temp, wei_temp;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        val_temp = Value[i];
        wei_temp = Weight[i];
        for(int j = wei_temp; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-wei_temp];
        }

        for(int j = 0; j < wei_temp; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        // 指定val为Value[i]
        val = _mm_set1_epi32(val_temp);
        for(int j = wei_temp; j <= Capacity; j+=4) {
            sol_cur = _mm_loadu_si128((__m128i*)& Ptr_sol_cur[j]);
            //val = _mm_loadu_si128((__m128i*)&Value[i]);
            sol_cur = _mm_add_epi32(sol_cur, val);

            sol_pre = _mm_loadu_si128((__m128i*)& Ptr_sol_pre[j]);
            sol_cur = _mm_max_epi32(sol_pre, sol_cur);
            _mm_storeu_si128((__m128i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "SSE Update Temp Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// SSE并行算法对齐内存
void Knapsack_SSE_update_temp_align(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m128i sol_pre, sol_cur, val;
    int val_temp, wei_temp;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        val_temp = Value[i];
        wei_temp = Weight[i];
        for(int j = wei_temp; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-wei_temp];
        }

        for(int j = 0; j < wei_temp; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        // 指定val为Value[i]
        val = _mm_set1_epi32(val_temp);
        if(wei_temp % 4 != 0){
            for(int j = wei_temp; j < wei_temp + 4 - wei_temp % 4; j++) {
                Ptr_sol_cur[j] += val_temp;
                Ptr_sol_cur[j] = (std::max)(Ptr_sol_pre[j], Ptr_sol_cur[j]);
            }
            wei_temp = wei_temp + 4 - wei_temp % 4;// 使得wei_temp为4的倍数，因为不再使用wei_temp所以可以直接修改
        }
        for(int j = wei_temp; j <= Capacity; j+=4) {
            sol_cur = _mm_load_si128((__m128i*)& Ptr_sol_cur[j]);
            //val = _mm_loadu_si128((__m128i*)&Value[i]);
            sol_cur = _mm_add_epi32(sol_cur, val);

            sol_pre = _mm_load_si128((__m128i*)& Ptr_sol_pre[j]);
            sol_cur = _mm_max_epi32(sol_pre, sol_cur);
            _mm_store_si128((__m128i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "SSE Update Temp Align Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// AVX并行算法
void Knapsack_AVX(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m256i sol_pre, sol_cur, val;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        for(int j = 0; j <= Capacity; j++) {
            if (j < Weight[i]) {
                Ptr_sol_cur[j] = -Value[i];
            }
            else {
                Ptr_sol_cur[j] = Ptr_sol_pre[j-Weight[i]];
            }
        }

        // 指定val为Value[i]
        val = _mm256_set1_epi32(Value[i]);
        for(int j = 0; j <= Capacity; j+=8) {
            sol_cur = _mm256_loadu_si256((__m256i*)& Ptr_sol_cur[j]);
            //val = _mm256_loadu_si256((__m256i*)&Value[i]);
            sol_cur = _mm256_add_epi32(sol_cur, val);

            sol_pre = _mm256_loadu_si256((__m256i*)& Ptr_sol_pre[j]);
            sol_cur = _mm256_max_epi32(sol_pre, sol_cur);
            _mm256_storeu_si256((__m256i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "AVX Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// AVX并行算法省略不需要的计算
void Knapsack_AVX_update(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m256i sol_pre, sol_cur, val;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        for(int j = Weight[i]; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-Weight[i]];
        }

        for(int j = 0; j < Weight[i]; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        // 指定val为Value[i]
        val = _mm256_set1_epi32(Value[i]);
        for(int j = Weight[i]; j <= Capacity; j+=8) {
            sol_cur = _mm256_loadu_si256((__m256i*)& Ptr_sol_cur[j]);
            //val = _mm256_loadu_si256((__m256i*)&Value[i]);
            sol_cur = _mm256_add_epi32(sol_cur, val);

            sol_pre = _mm256_loadu_si256((__m256i*)& Ptr_sol_pre[j]);
            sol_cur = _mm256_max_epi32(sol_pre, sol_cur);
            _mm256_storeu_si256((__m256i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "AVX Update Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// AVX并行算法提前设定临时变量存储Value[i]和Weight[i]避免运行中频繁访问数组
void Knapsack_AVX_update_temp(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m256i sol_pre, sol_cur, val;
    int val_temp, wei_temp;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        val_temp = Value[i];
        wei_temp = Weight[i];
        for(int j = wei_temp; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-wei_temp];
        }

        for(int j = 0; j < wei_temp; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        // 指定val为Value[i]
        val = _mm256_set1_epi32(val_temp);
        for(int j = wei_temp; j <= Capacity; j+=8) {
            sol_cur = _mm256_loadu_si256((__m256i*)& Ptr_sol_cur[j]);
            //val = _mm256_loadu_si256((__m256i*)&Value[i]);
            sol_cur = _mm256_add_epi32(sol_cur, val);

            sol_pre = _mm256_loadu_si256((__m256i*)& Ptr_sol_pre[j]);
            sol_cur = _mm256_max_epi32(sol_pre, sol_cur);
            _mm256_storeu_si256((__m256i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "AVX Update Temp Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// AVX并行算法对齐内存
void Knapsack_AVX_update_temp_align(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
    // 初始化
    int* Ptr_sol_pre = Solution_previous;
    int* Ptr_sol_cur = Solution_current;
    int* temp = nullptr;
    int Nsize = N;
    __m256i sol_pre, sol_cur, val;
    int val_temp, wei_temp;

    // 动态规划
    for(int i = 0; i < Nsize; i++) {
        val_temp = Value[i];
        wei_temp = Weight[i];
        for(int j = wei_temp; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
            Ptr_sol_cur[j] = Ptr_sol_pre[j-wei_temp];
        }

        for(int j = 0; j < wei_temp; j++){// 减少Weight[i]次计算和比较
            Ptr_sol_cur[j] = Ptr_sol_pre[j];
        }
        // 指定val为Value[i]
        val = _mm256_set1_epi32(val_temp);
        if(wei_temp % 8 != 0){
            for(int j = wei_temp; j < wei_temp + 8 - wei_temp % 8; j++) {
                Ptr_sol_cur[j] += val_temp;
                Ptr_sol_cur[j] = (std::max)(Ptr_sol_pre[j], Ptr_sol_cur[j]);
            }
            wei_temp = wei_temp + 8 - wei_temp % 8;// 使得wei_temp为8的倍数，因为不再使用wei_temp所以可以直接修改
        }
        for(int j = wei_temp; j <= Capacity; j+=8) {
            sol_cur = _mm256_load_si256((__m256i*)& Ptr_sol_cur[j]);
            //val = _mm256_loadu_si256((__m256i*)&Value[i]);
            sol_cur = _mm256_add_epi32(sol_cur, val);

            sol_pre = _mm256_load_si256((__m256i*)& Ptr_sol_pre[j]);
            sol_cur = _mm256_max_epi32(sol_pre, sol_cur);
            _mm256_store_si256((__m256i*)& Ptr_sol_cur[j], sol_cur);
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }

    // 输出结果
    //std::cout << "AVX Update Temp Align Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

int main()// Oo0LlIi
{
    // 内存计算 (2*Nsize*4 + 2*Capacity*4) B = 8*(Nsize + Capacity) B
    // 服务器 64KB 512KB 48MB
    // 8K 64K 6M
    // 本地 128KB(/2) 512KB(/2) 3MB
    // 16K(/2) 64K(/2) 384K
    int Nsize = 100;              // 数据规模
    int Capacity = 200000;          // 容量大小
    int Capacity_add_one = Capacity + 1;  // 数组大小
    int* Weight,* Value;          // 数据数组
    int* Solution_previous,* Solution_current; // 动态规划数组

    int repeat = 100;// 20;            // 重复次数 repeat=1,运行一个规模1分钟,总时间10分钟

    // 初始化计时器
    LARGE_INTEGER frequency, t1, t2;
    double elapsed_time = 0;
    QueryPerformanceFrequency(&frequency);// 获取计时器频率
    std::cout << "Frequency: " << frequency.QuadPart << std::endl;

    // 生成数据
    Weight = new int[Nsize];
    Value = new int[Nsize];
    Generate(Weight, Value, Nsize);
    Solution_current = (int*)_aligned_malloc(Capacity_add_one * sizeof(int), 32);
    Solution_previous = (int*)_aligned_malloc(Capacity_add_one * sizeof(int), 32);

    Zero(Solution_previous, Capacity_add_one);
    Zero(Solution_current, Capacity_add_one);

    for(int i = 0; i < repeat; i++){
        Zero(Solution_previous, Capacity_add_one);

        //Knapsack_Basic(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
        //Knapsack_Basic_update(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
        //Knapsack_Basic_update_temp(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
        Knapsack_SSE(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
        //Knapsack_SSE_update(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
        //Knapsack_SSE_update_temp(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
        //Knapsack_SSE_update_temp_align(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
        //其他
    }
    _aligned_free(Solution_current);
    _aligned_free(Solution_previous);
    delete[] Weight;
    delete[] Value;
    return 0;
}
