#include <bits/stdc++.h>
#include <omp.h>
#include <sys/time.h>
using namespace std;

// Structure for Item which stores weight and corresponding value of Item
struct Item {
    int weight;
    int value;
};

// Node structure to store information of decision tree
struct Node {
    int level, profit, bound;
    int weight;
};

#define W 600 // Weight of knapsack
#define n 100  // Number of items

Item items[n];
int maxProfit = 0;
omp_lock_t maxProfitLock;

// Comparison function to sort Item according to val/weight ratio
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

void knapsack_openmp(int tid, queue<Node>& Q) {
    while (!Q.empty()) {
        Node u = Q.front();
        Q.pop();

        if (u.level == -1)
            u.level = 0;

        if (u.level == n - 1)
            continue;

        Node v;
        v.level = u.level + 1;

        v.weight = u.weight + items[v.level].weight;
        v.profit = u.profit + items[v.level].value;

        if (v.weight <= W && v.profit > maxProfit) {
            omp_set_lock(&maxProfitLock);
            if (v.profit > maxProfit) {
                maxProfit = v.profit;
            }
            omp_unset_lock(&maxProfitLock);
        }

        v.bound = bound(v);

        if (v.bound > maxProfit)
            Q.push(v);

        v.weight = u.weight;
        v.profit = u.profit;
        v.bound = bound(v);
        if (v.bound > maxProfit)
            Q.push(v);
    }
}

int main() {
    srand(time(0));
    for (int i = 0; i < n; i++) {
        items[i].value = rand() % 100;
        items[i].weight = rand() % 100;
    }

    sort(items, items + n, cmp);

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

    queue<Node> threadQueues[4];
    for (int i = 0; i < 8; i++) {
        threadQueues[i % 4].push(nodes[i]);
    }

    omp_set_num_threads(4);
    omp_init_lock(&maxProfitLock);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        knapsack_openmp(tid, threadQueues[tid]);
    }

    gettimeofday(&end, NULL);
    double timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    timeuse /= 1000000;
    cout << "Maximum possible profit = " << maxProfit << endl;
    cout << "OpenMP time: " << timeuse << "s" << endl;

    maxProfit = 0;
    gettimeofday(&start, NULL);
    cout << "Maximum possible profit = " << knapsack() << endl;
    gettimeofday(&end, NULL);
    timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    timeuse /= 1000000;
    cout << "Serial time: " << timeuse << "s" << endl;

    omp_destroy_lock(&maxProfitLock);

    return 0;
}
