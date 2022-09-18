#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <atomic>
#include <chrono>
#include <iostream>
#include <vector>
#define THREADS 16
#define N 1000
using namespace std;
using namespace std::chrono;
const int M = 1000000;

int a[M];

int main() {
    int c;
    freopen("/Users/wayne_tx/Desktop/CG/out.txt", "w", stdout);
    for (c = 1; c <= 16; c++) {
        auto start1 = high_resolution_clock::now();
#pragma omp parallel for schedule(static) num_threads(c)
        for (int i = 0; i < M; i++) {
            a[i] = (8 * (i + 1) / 19 + 8) / c;
        }
        auto end1 = high_resolution_clock::now();
        auto ans1 = duration_cast<microseconds>(end1 - start1);
        printf("Number of Threads : %d\n", c);
        printf("Run Time : %d microseconds\n\n", ans1);
    }
    fclose(stdout);
    return 0;
}