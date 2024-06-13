#include <mpi.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector> /// 修改：使用优先队列
#include <queue>
// #include <priority_queue>
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
int localMaxProfit = 0; // 用来接收全局最大值
int myrank = 0, size = 0;
int one_finish = 0;
pthread_mutex_t maxProfitMutex = PTHREAD_MUTEX_INITIALIZER;

int allreduce_in_progress = 0; // 标记当前线程是否有未完成的allreduce

struct Node {
    int level = 0, profit = 0, bound = 0;
    int weight = 0;
    int lb = 0;/// 修改：增加表示当前节点的下界
};

class sortBylb {
public:
    bool operator()(const Node& a, const Node& b) const {/// 修改：按照下界从大到小排序
        return a.lb < b.lb;
    }
};

struct ThreadData { /// 修改：使用优先队列代替普通队列
    int threadId;
    // queue<Node> Q;
    priority_queue<Node, vector<Node>, sortBylb> Q;
};

bool cmp(Item a, Item b) {
    double r1 = (double)a.value / a.weight;
    double r2 = (double)b.value / b.weight;
    return r1 > r2;
}

void bound(Node& u) { /// 修改：不返回值，直接修改u的bound和lb
    if (u.weight > W)
        return; /// 去掉等于的情况

    int profit_bound = u.profit;
    int j = u.level;
    int totweight = u.weight;

    while ((j < n) && (totweight + items[j].weight <= W)) {
        totweight += items[j].weight;
        profit_bound += items[j].value;
        j++;
    }

    if (j < n) {/// 说明即将装入的物品不能完全装入背包
        u.lb = profit_bound;
        profit_bound += (W - totweight) * items[j].value / items[j].weight;
        u.bound = profit_bound;
    }
    else {
        u.lb = profit_bound;
        u.bound = profit_bound;
    }
    return;///
}

void *knapsack_pthread(void *threadData) {
    ThreadData *data = (ThreadData *)threadData;
    // queue<Node>& Q = data->Q;
    priority_queue<Node, vector<Node>, sortBylb>& Q = data->Q;
    Node u, v;
    int count = 0; /// 修改：由0改为-1，也即实现首次即同步，极大可能提前剪枝，不再修改因为已经初始化maxProfit
    MPI_Request request;

    // 计时
    struct timeval tstart_time, tend_time;
    gettimeofday(&tstart_time, NULL);

    while (!Q.empty()) {
        count++;
        u = Q.top();
        Q.pop();

        // 测试：输出节点u的信息
        cout << "rank " << myrank << " thread " << data->threadId << " level = " << u.level << " weight = " << u.weight << " profit = " << u.profit << " bound = " << u.bound << " lb = " << u.lb << endl;

        if (u.level >= n)
            continue;

        if (u.bound <= maxProfit) /// 重要 <=
            continue;/// 修改：提前剪枝

        /// 左节点
        v.weight = u.weight;
        v.profit = u.profit;
        v.level = u.level + 1; /// 发生错误的地方，未及时更新其level
        bound(v); /// 更新其上下界

        if (v.weight <= W && v.bound > maxProfit) {
            Q.push(v);
            if (v.lb > maxProfit) {
                pthread_mutex_lock(&maxProfitMutex);
                if (v.lb > maxProfit){
                    cout << "rank " << myrank << " thread " << data->threadId << " maxProfit change from " << maxProfit << " to " << v.lb << endl;
                    cout << "rank " << myrank << " thread " << data->threadId << " level = " << v.level << " weight = " << v.weight << " profit = " << v.profit << " bound = " << v.bound << " lb = " << v.lb << endl;
                    // maxProfit = v.bound;
                    maxProfit = v.lb;
                }
                pthread_mutex_unlock(&maxProfitMutex);
            }
        }


        /// 右节点
        v.level = u.level + 1;
        v.weight = u.weight + items[v.level - 1].weight; /// 修改：先计算右节点再计算左节点
        v.profit = u.profit + items[v.level - 1].value;
        bound(v); /// 更新其上下界

        if (v.weight <= W && v.bound > maxProfit) {
            Q.push(v);
            if (v.lb > maxProfit) {
                pthread_mutex_lock(&maxProfitMutex);
                if (v.lb > maxProfit){
                    cout << "rank " << myrank << " thread " << data->threadId << " maxProfit change from " << maxProfit << " to " << v.lb << endl;
                    cout << "rank " << myrank << " thread " << data->threadId << " level = " << v.level << " weight = " << v.weight << " profit = " << v.profit << " bound = " << v.bound << " lb = " << v.lb << endl;
                    // maxProfit = v.bound; /// 错误，应当是maxProfit = v.lb;
                    maxProfit = v.lb;
                }
                pthread_mutex_unlock(&maxProfitMutex);
            }
        }



        if(one_finish == 1){
            cout << "其他线程已经结束" << endl;
            /// continue;
        }

        if (data->threadId == 0 && size > 1) {
            // 发起新的allreduce
            if (!allreduce_in_progress && count % 100 == 0 && (maxProfit - localMaxProfit > 20)) { /// 1000 200
                cout << "rank " << myrank << " thread " << data->threadId << " allreduce start" << endl;
                MPI_Iallreduce(&maxProfit, &localMaxProfit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &request);
                allreduce_in_progress = 1;
            }

            // 检查未完成的allreduce
            if (allreduce_in_progress) {
                int flag;
                MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
                if (flag) {
                    allreduce_in_progress = 0;
                    cout << "R" << myrank << "T" << data->threadId << " allreduce finish!" << endl;
                    if (localMaxProfit > maxProfit){
                        cout << "Allreduce: R" << myrank << "T" << data->threadId << " maxProfit change from " << maxProfit << " to " << localMaxProfit << endl;
                        pthread_mutex_lock(&maxProfitMutex);
                        maxProfit = localMaxProfit;
                        pthread_mutex_unlock(&maxProfitMutex);
                    }
                }
            }
        }

        if (count % 10 == 0){ // 每10次检查一次线程状态
            cout << "TCheck: R" << myrank << "T" << data->threadId << " count = " << count << endl;
            cout << "TCheck: R" << myrank << "T" << data->threadId << " Q.size() = " << Q.size() << endl;
            cout << "TCheck: R" << myrank << "T" << data->threadId << " toplevel = " << Q.top().level << endl;
            cout << "TCheck: R" << myrank << "T" << data->threadId << " maxProfit = " << maxProfit << endl;
        }

    }

    gettimeofday(&tend_time, NULL);
    cout << "TFinal: R" << myrank << "T" << data->threadId << " calculation time: " << (tend_time.tv_sec - tstart_time.tv_sec) * 1000 + (tend_time.tv_usec - tstart_time.tv_usec) / 1000 << " ms" << endl;
    cout << "TFinal: R" << myrank << "T" << data->threadId << " maxProfit = " << maxProfit << endl;
    return NULL;
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
        // string filename = "data_ItemNum_100.txt";
        ifstream in(filename);
        if (!in) {
        cerr << "Unable to open file " << filename << endl;
        return 1;
        }
        for (int i = 0; i < n; i++) {
            in >> items[i].value >> items[i].weight;
        }
        in.close();
    }

    // 广播数据
    MPI_Bcast(items, n * sizeof(Item), MPI_BYTE, 0, MPI_COMM_WORLD);

    // knapsack_mpi
    gettimeofday(&start_time, NULL);
    priority_queue<Node, vector<Node>, sortBylb> Q;
    Node u, v;
    // 计算待分配节点的level和数量
    int nodelevel = int(ceil(log2(size)));
    int nodesize = int(pow(2, nodelevel));
    // 当前线程的节点数量和范围
    int localnodesize = int(ceil(nodesize / float(size)));
    int local_start = myrank * localnodesize;
    int local_end = min((myrank + 1) * localnodesize, nodesize);
    // 初始化进程节点
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
                // itemnow--; /// 错误发生的地方，应该在最后减一
            }
            itemnow--;
        }
        // if (nodelevel == 0) {
        //     u.weight = 0;
        //     u.profit = 0;
        // }

        bound(u); ///
        cout << "Rgenerate: R" << myrank << " level = " << u.level << " weight = " << u.weight << " profit = " << u.profit << " bound = " << u.bound << " lb = " << u.lb << endl;
        Q.push(u);
    }
    gettimeofday(&end_time, NULL);
    cout << "Rgenerate: R" << myrank << " node generation time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;
    cout << "Rgenerate: R" << myrank << " Q.size() = " << Q.size() << endl;
    gettimeofday(&start_time, NULL);

    // 创建线程节点
    pthread_t threads[NUM_THREADS];
    ThreadData threadData[NUM_THREADS];
    int itemsPerThread = int(ceil(Q.size() / float(NUM_THREADS)));
    // 根据线程数量生成子线程所需要的节点
    int rep = int(ceil(log2(itemsPerThread * NUM_THREADS / float(Q.size()))));
    int itemsAllThread = Q.size() * int(pow(2, rep));
    cout << "Tgenerate: R" << myrank << " itemsAllThread = " << itemsAllThread << " rep = " << rep << endl;
    // 确保Q中的节点均属于同一层
    for (int i = 0; i < rep; i++){
        int nowlevel = Q.top().level;
        Node temp;
        while (nowlevel == Q.top().level) {
            temp = Q.top();
            Q.pop();
            u.level = temp.level + 1;
            u.weight = temp.weight + items[u.level - 1].weight;
            u.profit = temp.profit + items[u.level - 1].value;
            bound(u); ///
            if (u.lb > maxProfit) maxProfit = u.lb; ///
            Q.push(u);

            v.level = temp.level + 1;
            v.weight = temp.weight;
            v.profit = temp.profit;
            bound(v); ///
            // if (u.lb > maxProfit) maxProfit = u.lb; /// 出现错误，应该是v.lb
            if (v.lb > maxProfit) maxProfit = v.lb; ///
            Q.push(v);
        }
    }
    cout << "Tgenerate: R" << myrank << " Q.size() = " << Q.size() << endl;


    itemsPerThread = int(floor(Q.size() / float(NUM_THREADS)));
    // 分配子线程的节点
    for (int i = 0; i < NUM_THREADS; i++) {
        if (i != NUM_THREADS - 1) {
            for (int j = 0; j < itemsPerThread; j++) {
                threadData[i].Q.push(Q.top());
                Q.pop();
                cout << "Tgenerate: R" << myrank << "T" << i << " bound " << threadData[i].Q.top().bound << " " << threadData[i].Q.top().lb << endl;
            }
        } else {
            while (!Q.empty()) {
                threadData[i].Q.push(Q.top());
                Q.pop();
                cout << "Tgenerate: R" << myrank << "T" << i << " bound " << threadData[i].Q.top().bound << " " << threadData[i].Q.top().lb << endl;
            }
        }
        threadData[i].threadId = i;
        // 输出每个子线程的节点数量
        cout << "Tgenerate: R" << myrank << "T" << i << " Q.size() = " << threadData[i].Q.size() << endl;
        pthread_create(&threads[i], NULL, knapsack_pthread, (void*)&threadData[i]);
    }

    // 等待线程结束
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    // if (Q.size() != 0) {
    //     cout << "Error: Q.size() != 0" << endl;
    // }

    // 广播本线程已经结束，用来停止同步
    // if(one_finish == 0){
    //     one_finish = 1;
    //     MPI_Bcast(&one_finish, 1, MPI_INT, myrank, MPI_COMM_WORLD);
    //     cout << "rank " << myrank << " one_finish = 1" << endl;
    // }

    gettimeofday(&end_time, NULL);
    cout << "RFinal: R" << myrank << " calculation time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;

    // 确保所有节点都计算完毕
    MPI_Barrier(MPI_COMM_WORLD);

    gettimeofday(&end_time, NULL);
    cout << "RFinal: R" << myrank << " total time: " << (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000 << " ms" << endl;

    // 计算全局最大值
    int globalMaxProfit = 0;
    MPI_Allreduce(&maxProfit, &globalMaxProfit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    cout <<"RFinal: R" << myrank << " localmaxprofit = " << maxProfit << endl;
    if (myrank == 0) {
        cout << "RFinal: globalMaxProfit = " << globalMaxProfit << endl;
    }

    MPI_Finalize();

    return 0;
}
