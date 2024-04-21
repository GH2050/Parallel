// 在Linux下编译运行，实现动态内存分配
// g++ -g -march=native ks_arm.cpp -o ks_arm
// 在Linux下进行模拟
// aarch64-linux-gnu-g++ -static -o neon -march=armv8.2-a neon.cpp 
// qemu−aarch64 neon

#include <iostream>

// Windows
//#include <windows.h> // 传统的min/max宏定义会影响程序中使用std::min/max，可以使用(std::min/max)阻止宏替换

// Linux
#include <time.h> // gettimeofday

// arm
#include <arm_neon.h> // NEON

// x86
//#include <nmmintrin.h> // SSE4.2
//#include <immintrin.h> // AVX2


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

// NEON并行算法
void Knapsack_NEON(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
   // 初始化
   int* Ptr_sol_pre = Solution_previous;
   int* Ptr_sol_cur = Solution_current;
   int* temp = nullptr;
   int Nsize = N;
   int32x4_t sol_pre, sol_cur, val;

   // 动态规划
   for(int i = 0; i < Nsize; i++) {
       for(int j = 0; j <= Capacity; j++) {
           if (j < Weight[i]) {
               Solution_current[j] = -Value[i];
           }
           else {
               Solution_current[j] = Solution_previous[j-Weight[i]];
           }
       }

       // 指定val为Value[i]
       val = vld1q_s32(&Value[i]);
       for(int j = 0; j <= Capacity; j+=4) {
           sol_cur = vld1q_s32(&Solution_current[j]);
           //val = vld1q_s32(&Value[i]);
           sol_cur = vaddq_s32(sol_cur, val);

           sol_pre = vld1q_s32(&Solution_previous[j]);
           sol_cur = vmaxq_s32(sol_pre, sol_cur);
           vst1q_s32(&Solution_current[j], sol_cur);
       }

       temp = Ptr_sol_pre;
       Ptr_sol_pre = Ptr_sol_cur;
       Ptr_sol_cur = temp;
   }

   // 输出结果
   //std::cout << "NEON Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// NEON并行算法省略不需要的计算
void Knapsack_NEON_update(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
   // 初始化
   int* Ptr_sol_pre = Solution_previous;
   int* Ptr_sol_cur = Solution_current;
   int* temp = nullptr;
   int Nsize = N;
   int32x4_t sol_pre, sol_cur, val;

   // 动态规划
   for(int i = 0; i < Nsize; i++) {
       for(int j = Weight[i]; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
           Solution_current[j] = Solution_previous[j-Weight[i]];
       }

       for(int j = 0; j < Weight[i]; j++){// 减少Weight[i]次计算和比较
           Ptr_sol_cur[j] = Ptr_sol_pre[j];
       }
       // 指定val为Value[i]
       val = vld1q_s32(&Value[i]);
       for(int j = Weight[i]; j <= Capacity; j+=4) {
           sol_cur = vld1q_s32(&Solution_current[j]);
           //val = vld1q_s32(&Value[i]);
           sol_cur = vaddq_s32(sol_cur, val);

           sol_pre = vld1q_s32(&Solution_previous[j]);
           sol_cur = vmaxq_s32(sol_pre, sol_cur);
           vst1q_s32(&Solution_current[j], sol_cur);
       }

       temp = Ptr_sol_pre;
       Ptr_sol_pre = Ptr_sol_cur;
       Ptr_sol_cur = temp;
   }

  // 输出结果
  //std::cout << "NEON Update Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// NEON并行算法提前设定临时变量存储Value[i]和Weight[i]避免运行中频繁访问数组
void Knapsack_NEON_update_temp(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
   // 初始化
   int* Ptr_sol_pre = Solution_previous;
   int* Ptr_sol_cur = Solution_current;
   int* temp = nullptr;
   int Nsize = N;
   int32x4_t sol_pre, sol_cur, val;
   int val_temp, wei_temp;

   // 动态规划
   for(int i = 0; i < Nsize; i++) {
       val_temp = Value[i];
       wei_temp = Weight[i];
       for(int j = wei_temp; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
           Solution_current[j] = Solution_previous[j-wei_temp];
       }

       for(int j = 0; j < wei_temp; j++){// 减少Weight[i]次计算和比较
           Ptr_sol_cur[j] = Ptr_sol_pre[j];
       }
       // 指定val为Value[i]
       val = vld1q_s32(&val_temp);
       for(int j = wei_temp; j <= Capacity; j+=4) {
           sol_cur = vld1q_s32(&Solution_current[j]);
           //val = vld1q_s32(&Value[i]);
           sol_cur = vaddq_s32(sol_cur, val);

           sol_pre = vld1q_s32(&Solution_previous[j]);
           sol_cur = vmaxq_s32(sol_pre, sol_cur);
           vst1q_s32(&Solution_current[j], sol_cur);
       }

       temp = Ptr_sol_pre;
       Ptr_sol_pre = Ptr_sol_cur;
       Ptr_sol_cur = temp;
   }

   // 输出结果
   //std::cout << "NEON Update Temp Result: " << Ptr_sol_pre[Capacity] << std::endl;
}

// NEON并行算法对齐内存
void Knapsack_NEON_update_temp_align(int* Weight, int* Value, int N, int Capacity, int* Solution_previous, int* Solution_current) {
   // 初始化
   int* Ptr_sol_pre = Solution_previous;
   int* Ptr_sol_cur = Solution_current;
   int* temp = nullptr;
   int Nsize = N;
   int32x4_t sol_pre, sol_cur, val;
   int val_temp, wei_temp;

   // 动态规划
   for(int i = 0; i < Nsize; i++) {
       val_temp = Value[i];
       wei_temp = Weight[i];
       for(int j = wei_temp; j <= Capacity; j++) {// 减少Weight[i]次计算赋值和Capability次判断
           Solution_current[j] = Solution_previous[j-wei_temp];
       }

       for(int j = 0; j < wei_temp; j++){// 减少Weight[i]次计算和比较
           Ptr_sol_cur[j] = Ptr_sol_pre[j];
       }
       // 指定val为Value[i]
       val = vld1q_s32(&val_temp);
       if(wei_temp % 4 != 0){
           for(int j = wei_temp; j < wei_temp + 4 - wei_temp % 4; j++) {
               Ptr_sol_cur[j] += val_temp;
               Ptr_sol_cur[j] = (std::max)(Ptr_sol_pre[j], Ptr_sol_cur[j]);
           }
           wei_temp = wei_temp + 4 - wei_temp % 4;// 使得wei_temp为4的倍数，因为不再使用wei_temp所以可以直接修改
       }
       for(int j = wei_temp; j <= Capacity; j+=4) {
           sol_cur = vld1q_s32(&Solution_current[j]);
           //val = vld1q_s32(&Value[i]);
           sol_cur = vaddq_s32(sol_cur, val);

           sol_pre = vld1q_s32(&Solution_previous[j]);
           sol_cur = vmaxq_s32(sol_pre, sol_cur);
           vst1q_s32(&Solution_current[j], sol_cur);
       }

       temp = Ptr_sol_pre;
       Ptr_sol_pre = Ptr_sol_cur;
       Ptr_sol_cur = temp;
   }

   // 输出结果
   //std::cout << "NEON Update Temp Align Result: " << Ptr_sol_pre[Capacity] << std::endl;
}
   
int main()// Oo0LlIi
{
   // 内存计算 (2*Nsize*4 + 2*Capacity*4) B = 8*(Nsize + Capacity) B
   // 服务器 64KB 512KB 48MB
   // 8K 64K 6M
   // 本地 128KB(/2) 512KB(/2) 3MB
   // 16K(/2) 64K(/2) 384K
   int Nsize = 100;              // 数据规模
   int Capacity = 0;          // 容量大小
   int Capacity_add_one = Capacity + 1;  // 数组大小
   int* Weight,* Value;          // 数据数组
   int* Solution_previous,* Solution_current; // 动态规划数组
   
   int repeat = 40;// 20;            // 重复次数
   
   // 初始化计时器
   struct timespec t1, t2;       // 计时器
   time_t elapsed_time = 0;       // 时间，单位为ms
   
   // 生成数据
   Weight = new int[Nsize];
   Value = new int[Nsize];
   Generate(Weight, Value, Nsize);

    //for (Capacity = 500000; Capacity <= 8000000; Capacity += 500000) {
    //for (Capacity = 100000; Capacity <= 500000; Capacity += 50000) {
    for(Capacity = 2000; Capacity <= 100000; Capacity += 2000){
       std::cout << "Capacity:" << Capacity << std::endl;
       Capacity_add_one = Capacity + 1;
       Solution_current = (int *)aligned_alloc(32, Capacity_add_one * sizeof(int));
       Solution_previous = (int *)aligned_alloc(32, Capacity_add_one * sizeof(int));
       Zero(Solution_previous, Capacity_add_one);
       Zero(Solution_current, Capacity_add_one);

       // 串行算法
       for(int i = 0; i < repeat; i++){
           Zero(Solution_previous, Capacity_add_one);
           clock_gettime(CLOCK_REALTIME, &t1);
           Knapsack_Basic(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
           clock_gettime(CLOCK_REALTIME, &t2);
           elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000;
       }
       std::cout << "Basic Time: " << elapsed_time / 1000.0 / repeat << "ms" << std::endl;
       elapsed_time = 0;

       // 串行算法省略不需要的计算
       for(int i = 0; i < repeat; i++){
           Zero(Solution_previous, Capacity_add_one);
           clock_gettime(CLOCK_REALTIME, &t1);
           Knapsack_Basic_update(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
           clock_gettime(CLOCK_REALTIME, &t2);
           elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000;
       }
       std::cout << "Basic Update Time: " << elapsed_time / 1000.0 / repeat << "ms" << std::endl;
       elapsed_time = 0;

       // 串行算法提前设定临时变量存储Value[i]和Weight[i]避免运行中频繁访问数组
       for(int i = 0; i < repeat; i++){
           Zero(Solution_previous, Capacity_add_one);
           clock_gettime(CLOCK_REALTIME, &t1);
           Knapsack_Basic_update_temp(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
           clock_gettime(CLOCK_REALTIME, &t2);
           elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000;
       }
       std::cout << "Basic Update Temp Time: " << elapsed_time / 1000.0 / repeat << "ms" << std::endl;
       elapsed_time = 0;

       // NEON
       for(int i = 0; i < repeat; i++){
           Zero(Solution_previous, Capacity_add_one);
           clock_gettime(CLOCK_REALTIME, &t1);
           Knapsack_NEON(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
           clock_gettime(CLOCK_REALTIME, &t2);
           elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000;
       }
       std::cout << "NEON Time: " << elapsed_time / 1000.0 / repeat << "ms" << std::endl;
       elapsed_time = 0;

       // NEON省略不需要的计算
       for(int i = 0; i < repeat; i++){
           Zero(Solution_previous, Capacity_add_one);
           clock_gettime(CLOCK_REALTIME, &t1);
           Knapsack_NEON_update(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
           clock_gettime(CLOCK_REALTIME, &t2);
           elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000;
       }
       std::cout << "NEON Update Time: " << elapsed_time / 1000.0 / repeat << "ms" << std::endl;
       elapsed_time = 0;

       // NEON提前设定临时变量存储Value[i]和Weight[i]避免运行中频繁访问数组
       for(int i = 0; i < repeat; i++){
           Zero(Solution_previous, Capacity_add_one);
           clock_gettime(CLOCK_REALTIME, &t1);
           Knapsack_NEON_update_temp(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
           clock_gettime(CLOCK_REALTIME, &t2);
           elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000;
       }
       std::cout << "NEON Update Temp Time: " << elapsed_time / 1000.0 / repeat << "ms" << std::endl;
       elapsed_time = 0;

       // NEON对齐内存
       for(int i = 0; i < repeat; i++){
           Zero(Solution_previous, Capacity_add_one);
           clock_gettime(CLOCK_REALTIME, &t1);
           Knapsack_NEON_update_temp_align(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);
           clock_gettime(CLOCK_REALTIME, &t2);
           elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_nsec - t1.tv_nsec) / 1000;
       }
       std::cout << "NEON Update Temp Align Time: " << elapsed_time / 1000.0 / repeat << "ms" << std::endl;
       elapsed_time = 0;
       
		//aligned_free(Solution_current);
		//aligned_free(Solution_previous);
       free(Solution_current);
       free(Solution_previous);
   }
   delete[] Weight;
   delete[] Value;
   return 0;
}