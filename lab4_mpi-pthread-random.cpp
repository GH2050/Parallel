#include <mpi.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <queue>
#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <string>
#include <pthread.h>
using namespace std;

#define W 800 // Weight of knapsack
#define n 100 // Number of items
#define NUM_THREADS 2 // Number of threads

struct Item {
    int weight = 0;
    int value = 0;
};

Item items[n];
int maxProfit = 0;
int localMaxProfit = 0;// 用来接收全局最大值
int myrank = 0, size = 0;
int one_finish = 0;
pthread_mutex_t maxProfitMutex = PTHREAD_MUTEX_INITIALIZER;

int allreduce_in_progress = 0; // 标记当前线程是否有未完成的allreduce

struct Node {
    int level = 0, profit = 0, bound = 0;
    int weight = 0;
};

struct ThreadData {
    int threadId;
    queue<Node> Q;
};

bool cmp(Item a, Item b) {
    double r1 = (double)a.value / a.weight;
    double r2 = (double)b.value / b.weight;
    return r1 > r2;
}

int bound(Node& u) { // 需要修改
    if (u.weight >= W)
        return 0;

    int profit_bound = u.profit;
    int j = u.level;
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

void *knapsack_pthread(void *threadData) {
    ThreadData *data = (ThreadData *)threadData;
    queue<Node>& Q = data->Q;
    Node u, v;
    int count = 0;
    MPI_Request request;

    // 计时
    struct timeval tstart_time, tend_time;
    gettimeofday(&tstart_time, NULL);

    while (!Q.empty()) {
        count++;
        u = Q.front();
        Q.pop();

        if (u.level == -1)
            v.level = 0;

        if (u.level >= n)
            continue;

        v.level = u.level + 1;
        v.weight = u.weight + items[v.level - 1].weight;
        v.profit = u.profit + items[v.level - 1].value;

        if (v.weight <= W && v.profit > maxProfit) {
            pthread_mutex_lock(&maxProfitMutex);
            if (v.profit > maxProfit){
                maxProfit = v.profit;
                // cout << "rank " << myrank << " thread " << data->threadId << " maxProfit change to " << maxProfit << endl;
            }
            pthread_mutex_unlock(&maxProfitMutex);
        }

        // if (v.weight <= W && v.profit > maxProfit) {
        //     pthread_mutex_lock(&maxProfitMutex);
        //     if (v.profit > maxProfit)
        //         maxProfit = v.profit;
        //     pthread_mutex_unlock(&maxProfitMutex);
        // }
        if(int t = bound(v)){
            v.bound = t;
            if (v.bound > maxProfit){
                Q.push(v);
                // cout << "Push: " << v.bound << ">" << maxProfit << endl;
            }
        }



        v.weight = u.weight;
        v.profit = u.profit;
        if(int t = bound(v))
            v.bound = t;
        if (v.bound > maxProfit)
            Q.push(v);

        if(one_finish == 1){
            cout << "其他线程已经结束" << endl;
            continue;
        }

        if (data->threadId == 0 && size > 1) {
            // 非阻塞的Allreduce
            if (count % 1000 == 0 && !allreduce_in_progress && (maxProfit - localMaxProfit > 200)) {
                cout << "rank " << myrank << " thread " << data->threadId << " allreduce" << endl;
                MPI_Iallreduce(&maxProfit, &localMaxProfit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &request);
                allreduce_in_progress = 1;
            }

            // 检查是否有未完成的allreduce
            if (allreduce_in_progress) {
                int flag;
                MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
                if (flag) {
                    allreduce_in_progress = 0;
                    cout << "rank " << myrank << " thread " << data->threadId << " allreduce finish" << endl;
                    if (localMaxProfit > maxProfit){
                        cout << "allreduce maxProfit change!" << endl;
                        cout << "rank " << myrank << " thread " << data->threadId << " maxProfit change to " << localMaxProfit << endl;
                        pthread_mutex_lock(&maxProfitMutex);
                        maxProfit = localMaxProfit;
                        pthread_mutex_unlock(&maxProfitMutex);
                    }
                }
            }
        }

        if (count % 1000000 == 0){
            cout << "rank " << myrank << " thread " << data->threadId << " count = " << count << endl;
            cout << Q.size() << endl;
            cout << Q.front().level << endl;
            cout << "maxP = " << maxProfit << endl;
        }

    }

    gettimeofday(&tend_time, NULL);
    cout << "rank " << myrank <<" Thread " << data->threadId << " calculation time: " << (tend_time.tv_sec - tstart_time.tv_sec) * 1000 + (tend_time.tv_usec - tstart_time.tv_usec) / 1000 << " ms" << endl;

    return NULL;
}








int main(int argc, char* argv[]) {
    struct timeval start_time, end_time;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    cout<<"rank "<<myrank<<" size "<<size<<endl;

    // 生成数据
    if (myrank == 0) {
        gettimeofday(&start_time, NULL);
        for (int i = 0; i < n; i++) {
            items[i].value = rand() % 100 + 1;
            items[i].weight = rand() % 100 + 1;
        }
        sort(items, items + n, cmp);
        // 保存该数据以便下次读取，命名形式为data_ItemNum_n.txt
        string filename = "data_ItemNum_" + to_string(n) + ".txt";
        ofstream out(filename);
        for (int i = 0; i < n; i++) {
            out << items[i].value << " " << items[i].weight << endl;
        }
        out.close();
        gettimeofday(&end_time, NULL);
        cout << "Data generation time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;
    }

    // // 读取数据
    // if (myrank == 0) {
    //     string filename = "data_ItemNum_" + to_string(n) + ".txt";
    //     // string filename = "data_ItemNum_100.txt";
    //     ifstream in(filename);
    //     for (int i = 0; i < n; i++) {
    //         in >> items[i].value >> items[i].weight;
    //     }
    //     in.close();
    // }

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
        if (nodelevel == 0) {
            u.weight = 0;
            u.profit = 0;
        }

        if(int t =bound(u))
            u.bound = t;
        cout << "rank " << myrank << " u.level = " << u.level << " u.weight = " << u.weight << " u.profit = " << u.profit << " u.bound = " << u.bound << endl;
        Q.push(u);
    }
    gettimeofday(&end_time, NULL);
    cout << "rank " << myrank << " node generation time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;
    cout << "rank " << myrank << " Q.size() = " << Q.size() << endl;
    gettimeofday(&start_time, NULL);

    // 创建线程
    pthread_t threads[NUM_THREADS];
    ThreadData threadData[NUM_THREADS];
    int itemsPerThread = int(ceil(Q.size() / float(NUM_THREADS)));
    // 根据线程数量生成子线程所需要的节点
    int rep = int(ceil(log2(itemsPerThread * NUM_THREADS / float(Q.size()))));
    int itemsAllThread = Q.size() * pow(2, rep);
    cout << "rank " << myrank << " itemsAllThread = " << itemsAllThread << " rep = " << rep << endl;
    // 确保Q中的节点均属于同一层
    for (int i = 0; i < rep; i++){
        int nowlevel = Q.front().level;
        Node temp;
        while (nowlevel == Q.front().level) {
            temp = Q.front();
            Q.pop();
            u.level = temp.level + 1;
            u.weight = temp.weight + items[u.level - 1].weight;
            u.profit = temp.profit + items[u.level - 1].value;
            if(int t =bound(v))
                v.bound = t;
            Q.push(u);
            v.level = temp.level + 1;
            v.weight = temp.weight;
            v.profit = temp.profit;
            v.bound = temp.bound;
            Q.push(v);
        }
    }
    cout << "rank " << myrank << " Q.size() = " << Q.size() << endl;
    itemsPerThread = int(floor(Q.size() / float(NUM_THREADS)));
    // 生成子线程的节点
    for (int i = 0; i < NUM_THREADS; i++) {
        if (i != NUM_THREADS - 1) {
            for (int j = 0; j < itemsPerThread; j++) {
                threadData[i].Q.push(Q.front());
                Q.pop();
            }
        } else {
            while (!Q.empty()) {
                threadData[i].Q.push(Q.front());
                Q.pop();
            }
        }
        threadData[i].threadId = i;
        // 输出每个子线程的节点数量
        cout << "rank " << myrank << " thread " << i << " Q.size() = " << threadData[i].Q.size() << endl;
        pthread_create(&threads[i], NULL, knapsack_pthread, (void*)&threadData[i]);
    }

    // 等待线程结束
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    if (Q.size() != 0) {
        cout << "Error: Q.size() != 0" << endl;
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
