#include <mpi.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <queue>
#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <string>
using namespace std;

#define W 1000 // Weight of knapsack
#define n 100 // Number of items

struct Item {
    int weight;
    int value;
};

Item items[n];
int maxProfit = 0;
int localMaxProfit = 0;// 用来接收全局最大值
int myrank = 0, size = 0;
int one_finish = 0;

struct Node {
    int level, profit, bound;
    int weight;
};

bool cmp(Item a, Item b) {
    double r1 = (double)a.value / a.weight;
    double r2 = (double)b.value / b.weight;
    return r1 > r2;
}

int bound(Node& u) {
    if (u.weight >= W)
        return 0;

    int profit_bound = u.profit;
    int j = u.level + 1;
    int totweight = u.weight;

    while ((j < n) && (totweight + items[j].weight <= W)) {
        totweight += items[j].weight;
        profit_bound += items[j].value;
        j++;
    }

    if (j < n)
        profit_bound += (W - totweight) * items[j].value / items[j].weight;

    return profit_bound;
}

int main(int argc, char* argv[]) {
    struct timeval start_time, end_time;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    cout<<"rank "<<myrank<<" size "<<size<<endl;

    // 生成数据    
    // if (myrank == 0) {
    //     gettimeofday(&start_time, NULL);
    //     for (int i = 0; i < n; i++) {
    //         items[i].value = rand() % 100 + 1;
    //         items[i].weight = rand() % 100 + 1;
    //     }
    //     sort(items, items + n, cmp);
    //     // 保存该数据以便下次读取，命名形式为data_ItemNum_n.txt
    //     string filename = "data_ItemNum_" + to_string(n) + ".txt";
    //     ofstream out(filename);
    //     for (int i = 0; i < n; i++) {
    //         out << items[i].value << " " << items[i].weight << endl;
    //     }
    //     out.close();
    //     gettimeofday(&end_time, NULL);
    //     cout << "Data generation time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;
    // }

    // 读取数据
    if (myrank == 0) {
        string filename = "data_ItemNum_" + to_string(n) + ".txt";
        ifstream in(filename);
        for (int i = 0; i < n; i++) {
            in >> items[i].value >> items[i].weight;
        }
        in.close();
    }
    
    // 广播数据
    MPI_Bcast(items, n * sizeof(Item), MPI_BYTE, 0, MPI_COMM_WORLD);

    // knapsack_mpi
    gettimeofday(&start_time, NULL);
    queue<Node> Q;
    Node u, v;
    // 计算待分配节点的level和数量
    int nodelevel = int(ceil(log2(size)));
    int nodesize = int(pow(2, nodelevel));
    // 当前线程的节点数量和范围
    int localnodesize = int(ceil(nodesize / float(size)));
    int local_start = myrank * localnodesize;
    int local_end = min((myrank + 1) * localnodesize, nodesize);
    // 初始化节点
    for (int i = local_start; i < local_end; i++) {
        u.level = nodelevel;
        int quotient = i;
        int remainder = 0;
        int itemnow = nodelevel - 1;
        while (quotient > 0) {
            remainder = quotient % 2;
            quotient = quotient / 2;
            if (remainder == 1) {
                u.weight += items[itemnow].weight;
                u.profit += items[itemnow].value;
                itemnow--;
            }
        }
        u.bound = bound(u);
        Q.push(u);
    }
    gettimeofday(&end_time, NULL);
    cout << "rank " << myrank << " node generation time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;
    gettimeofday(&start_time, NULL);
    int count = 0;
    MPI_Request request;
    int allreduce_in_progress = 0; // 标记当前线程是否有未完成的allreduce
    // 进行计算
    while (!Q.empty()) {
        count++;
        u = Q.front();
        Q.pop();

        if (u.level == -1)
            v.level = 0;

        if (u.level == n - 1)
            continue;

        v.level = u.level + 1;

        v.weight = u.weight + items[v.level].weight;
        v.profit = u.profit + items[v.level].value;

        if (v.weight <= W && v.profit > maxProfit) {
            maxProfit = v.profit;
        }

        v.bound = bound(v);

        if (v.bound > maxProfit)
            Q.push(v);

        v.weight = u.weight;
        v.profit = u.profit;
        v.bound = bound(v);
        if (v.bound > maxProfit)
            Q.push(v);

        if(one_finish == 1){
            cout << "其他线程已经结束" << endl;
            continue;
        }

        // 非阻塞的Allreduce
        if (count % 2000 == 0 && !allreduce_in_progress && (maxProfit - localMaxProfit > 200)) {
            cout << "rank " << myrank << " allreduce" << endl;
            MPI_Iallreduce(&maxProfit, &localMaxProfit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &request);
            allreduce_in_progress = 1;
        }

        // 检查是否有未完成的allreduce
        if (allreduce_in_progress) {
            int flag;
            MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
            if (flag) {
                allreduce_in_progress = 0;
                cout << "rank " << myrank << " allreduce finish" << endl;
                if (localMaxProfit > maxProfit){
                    cout << "rank " << myrank << " maxProfit change to " << localMaxProfit << endl;
                    maxProfit = localMaxProfit;
                }
            }
        }
    }

    // 广播本线程已经结束，用来停止同步
    if(one_finish == 0){
        one_finish = 1;
        MPI_Bcast(&one_finish, 1, MPI_INT, myrank, MPI_COMM_WORLD);
        cout << "rank " << myrank << " one_finish = 1" << endl;
    }

    gettimeofday(&end_time, NULL);
    cout << "rank " << myrank << " calculation time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;
    
    // 确保所有节点都计算完毕
    MPI_Barrier(MPI_COMM_WORLD);
    
    gettimeofday(&end_time, NULL);
    cout << "rank " << myrank << " total time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;

    // 计算全局最大值
    int globalMaxProfit = 0;
    MPI_Allreduce(&maxProfit, &globalMaxProfit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (myrank == 0) {
        cout << "Maximum possible profit = " << globalMaxProfit << endl;
    }

    MPI_Finalize();

    return 0;
}
