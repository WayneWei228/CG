#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <vector>
#include <atomic>
#include <mutex>
#define THREADS 16
#define N 1000
#define CHUNK 100
using namespace std;

atomic_int total(0);

struct Coordinate {
    int X, Y;

    bool operator==(const Coordinate& that) const { return X == that.X && Y == that.Y; }
};

struct CoordinateHash {
    size_t operator()(const Coordinate& that) const { return that.X * size_t(9875321) + that.Y; }
};


int main() {
    clock_t start, end;
    start = clock();
    vector<int> a;
    int i;
    mutex mx;
#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
    for (i = 0; i < N; i++) {
        mx.lock();
        a.emplace_back(i);
        mx.unlock();
    }
    printf("%d\n", int(a.size()));
    /* all threads done */
    printf("All done!\n");
    printf("Total %u\n", total.load());
    end = clock();
    return 0;
}