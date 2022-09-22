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
    int c = 8;

    // for (c = 1; c <= 16; c++) {
        printf("Number of Threads : %d\n", c);
        sleep(3);
        auto start1 = steady_clock::now();
#pragma omp parallel for schedule(static) num_threads(c)
        for (long int i = 0; i < 10000000000; i++) {
            a[i % M] = i;
        }
        auto end1 = steady_clock::now();
        auto ans1 = duration_cast<milliseconds>(end1 - start1).count();
        printf("Run Time : %lld milliseconds\n\n", ans1);
    return 0;
}