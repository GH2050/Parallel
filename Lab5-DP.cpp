#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

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
                Ptr_sol_cur[j] = Ptr_sol_pre[j];
            } else {
                Ptr_sol_cur[j] = std::max(Ptr_sol_pre[j], Ptr_sol_pre[j - Weight[i]] + Value[i]);
            }
        }

        temp = Ptr_sol_pre;
        Ptr_sol_pre = Ptr_sol_cur;
        Ptr_sol_cur = temp;
    }
}


int main() {
    int Nsize = 100;              // 数据规模
    int Capacity = 800;           // 容量大小

    int Weight[100]{}, Value[100]{}; // 数据数组
    int Solution_previous[801]{}, Solution_current[801]{}; // 动态规划数组

    // 读取数据
    string filename = "data_ItemNum_100.txt";
    ifstream in(filename);
    if (!in) {
        cerr << "Unable to open file " << filename << endl;
        return 1;
    }
    for (int i = 0; i < Nsize; i++) {
        in >> Value[i] >> Weight[i];
    }
    in.close();


    // 动态规划算法
    Knapsack_Basic(Weight, Value, Nsize, Capacity, Solution_previous, Solution_current);

    // 输出结果
    cout << "Final Result: " << Solution_previous[Capacity] << endl;

    return 0;
}
