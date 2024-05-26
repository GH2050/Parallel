#include <bits/stdc++.h>
#include <algorithm>
#include <queue>
#include <iostream>
using namespace std;

// Structure for Item which store weight and corresponding
// value of Item
struct Item
{
	int weight;
	int value;
};

// Node structure to store information of decision
// tree
struct Node
{
	// level --> Level of node in decision tree (or index
	//			 in arr[]
	// profit --> Profit of nodes on path from root to this
	//		 node (including this node)
	// bound ---> Upper bound of maximum profit in subtree
	//		 of this node/
	int level, profit, bound;
	int weight;
};

// 使用pthread实现静态分配，每个子线程包含一个队列
#include <pthread.h>
#include <sys/time.h>

#define NUM_THREADS 6 // Number of threads
#define W 800 // Weight of knapsack // Capacity of knapsack 2000
#define n 100 // Number of items

struct ThreadData
{
	int threadId;
	queue<Node> Q;
};

Item items[n];
int maxProfit = 0;
pthread_mutex_t maxProfitMutex = PTHREAD_MUTEX_INITIALIZER;

// Comparison function to sort Item according to
// val/weight ratio
bool cmp(Item a, Item b)
{
	double r1 = (double)a.value / a.weight;
	double r2 = (double)b.value / b.weight;
	return r1 > r2;
}

int bound(Node& u)
{
	// if weight overcomes the knapsack capacity, return
	// 0 as expected bound
	if (u.weight >= W)
		return 0;

	// initialize bound on profit by current profit
	int profit_bound = u.profit;

	// start including items from index 1 more to current
	// item index
	int j = u.level + 1;
	int totweight = u.weight;

	// checking index condition and knapsack capacity
	// condition
	while ((j < n) && (totweight + items[j].weight <= W))
	{
		totweight += items[j].weight;
		profit_bound += items[j].value;
		j++;
	}

	// If j is not n, include last item partially for
	// upper bound on profit
	if (j < n)
		profit_bound += (W - totweight) * items[j].value /
		items[j].weight;

	return profit_bound;
}

void *knapsack_pthread(void *threadData)
{
	ThreadData *data = (ThreadData *)threadData;
	int tid = data->threadId;
	queue<Node>& Q = data->Q;

	// 增加计时功能
	struct timeval start, end;
	gettimeofday(&start, NULL);

	// One by one extract an item from decision tree
	// compute profit of all children of extracted item
	// and keep saving maxProfit
	while (!Q.empty())
	{
		// Dequeue a node
		Node u = Q.front();
		Q.pop();

		// If it is starting node, assign level 0
		if (u.level == -1)
			u.level = 0;

		// If there is nothing on next level
		if (u.level == n - 1)
			continue;

		// Else if not last node, then increment level,
		// and compute profit of children nodes.
		Node v;
		v.level = u.level + 1;

		// Taking current level's item add current
		// level's weight and value to node u's
		// weight and value
		v.weight = u.weight + items[v.level].weight;
		v.profit = u.profit + items[v.level].value;

		// If cumulated weight is less than W and
		// profit is greater than previous profit,
		// update maxprofit
		if (v.weight <= W && v.profit > maxProfit)
		{
			pthread_mutex_lock(&maxProfitMutex);
			if (v.profit > maxProfit)
			{
				maxProfit = v.profit;
			}
			pthread_mutex_unlock(&maxProfitMutex);
		}

		// Get the upper bound on profit to decide
		// whether to add v to Q or not.
		v.bound = bound(v);// n, W, items均为全局变量

		// If bound value is greater than profit,
		// then only push into queue for further
		// consideration
		if (v.bound > maxProfit)
			Q.push(v);

		// Do the same thing, but Without taking the item in knapsack
		v.weight = u.weight;
		v.profit = u.profit;
		v.bound = bound(v);
		if (v.bound > maxProfit)
			Q.push(v);
	}

	// 计时结束
	gettimeofday(&end, NULL);
	double timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	timeuse /= 1000000;
	std::cout << "Thread " << tid << " time: " << timeuse << "s" << std::endl;

	return NULL;
}

// 串行算法
int knapsack(){
	// make a queue for traversing the node
	queue<Node> Q;
	Node u, v;

	// dummy node at starting
	u.level = -1;
	u.profit = u.weight = 0;
	Q.push(u);

	// One by one extract an item from decision tree
	// compute profit of all children of extracted item
	// and keep saving maxProfit
	while (!Q.empty())
	{
		// Dequeue a node
		u = Q.front();
		Q.pop();

		// If it is starting node, assign level 0
		if (u.level == -1)
			v.level = 0;

		// If there is nothing on next level
		if (u.level == n - 1)
			continue;

		// Else if not last node, then increment level,
		// and compute profit of children nodes.
		v.level = u.level + 1;

		// Taking current level's item add current
		// level's weight and value to node u's
		// weight and value
		v.weight = u.weight + items[v.level].weight;
		v.profit = u.profit + items[v.level].value;

		// If cumulated weight is less than W and
		// profit is greater than previous profit,
		// update maxprofit
		if (v.weight <= W && v.profit > maxProfit)
			maxProfit = v.profit;

		// Get the upper bound on profit to decide
		// whether to add v to Q or not.
		v.bound = bound(v);// n, W, items均为全局变量

		// If bound value is greater than profit,
		// then only push into queue for further
		// consideration
		if (v.bound > maxProfit)
			Q.push(v);

		// Do the same thing, but Without taking the item in knapsack
		v.weight = u.weight;
		v.profit = u.profit;
		v.bound = bound(v);
		if (v.bound > maxProfit)
			Q.push(v);
	}

	return maxProfit;

}

int main()
{

	// 随机生成物品的价值和重量
	for (int i = 0; i < n; i++)
	{
		items[i].value = rand() % 100;
		items[i].weight = rand() % 100;
	}

	// sorting Item on basis of value per unit weight.
	sort(items, items + n, cmp);

	// 创建8个Node用于初始化子线程的队列，将前三个物品是否放入背包分成8个子问题
	Node nodes[8];
	nodes[0].level = 2;
	nodes[0].profit = items[0].value + items[1].value + items[2].value;
	nodes[0].weight = items[0].weight + items[1].weight + items[2].weight;
	nodes[0].bound = bound(nodes[0]);
	nodes[1].level = 2;
	nodes[1].profit = items[0].value + items[1].value;
	nodes[1].weight = items[0].weight + items[1].weight;
	nodes[1].bound = bound(nodes[1]);
	nodes[2].level = 2;
	nodes[2].profit = items[0].value + items[2].value;
	nodes[2].weight = items[0].weight + items[2].weight;
	nodes[2].bound = bound(nodes[2]);
	nodes[3].level = 2;
	nodes[3].profit = items[1].value + items[2].value;
	nodes[3].weight = items[1].weight + items[2].weight;
	nodes[3].bound = bound(nodes[3]);
	nodes[4].level = 2;
	nodes[4].profit = items[0].value;
	nodes[4].weight = items[0].weight;
	nodes[4].bound = bound(nodes[4]);
	nodes[5].level = 2;
	nodes[5].profit = items[1].value;
	nodes[5].weight = items[1].weight;
	nodes[5].bound = bound(nodes[5]);
	nodes[6].level = 2;
	nodes[6].profit = items[2].value;
	nodes[6].weight = items[2].weight;
	nodes[6].bound = bound(nodes[6]);
	nodes[7].level = 2;
	nodes[7].profit = 0;
	nodes[7].weight = 0;
	nodes[7].bound = bound(nodes[7]);

	// 创建线程数据
	ThreadData threadData[NUM_THREADS];
	for (int i = 0; i < NUM_THREADS; i++)
	{
		threadData[i].threadId = i;
	}
	for (int i = 0; i < 8; i++)
	{
		threadData[i % NUM_THREADS].Q.push(nodes[i]);
	}

	std::cout << "Main thread started" << endl;

	// 增加计时功能
	struct timeval start, end;
	gettimeofday(&start, NULL);

	// 创建线程
	pthread_t threads[NUM_THREADS];
	int rc;
	int i;
	for (i = 0; i < NUM_THREADS; i++)
	{
		rc = pthread_create(&threads[i], NULL, knapsack_pthread, (void *)&threadData[i]);
		if (rc)
		{
			std::cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}

	// 等待线程结束
	for (i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// 计时结束
	gettimeofday(&end, NULL);
	double timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	timeuse /= 1000000;
	std::cout << "Maximum possible profit = " << maxProfit << endl;
	std::cout << "Main thread time: " << timeuse << "s" << std::endl;

	// 作为对比，串行计算
	// maxProfit = 0;
	// gettimeofday(&start, NULL);
	// std::cout << "Maximum possible profit = " << knapsack() << endl;
	// gettimeofday(&end, NULL);
	// timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	// timeuse /= 1000000;
	// std::cout << "Serial time: " << timeuse << "s" << std::endl;

	return 0;
}
