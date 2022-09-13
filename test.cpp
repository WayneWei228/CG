#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <iostream>
#include <mutex>
#include <vector>
#define THREADS 16
#define N 1000
#define CHUNK 100
using namespace std;
using namespace std::chrono;
const int M = 1000000;

int main() {
    vector<int> A;
    A.resize(M);
    for (int num = 1; num <= 16; num++) {
        auto start = high_resolution_clock::now();
#pragma omp parallel for schedule(static) num_threads(num)
        for (int i = 0; i < A.size(); i++) {
            A[i] = i * 2 - 1;
        }

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Number of threads : " << num << endl;
        cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
        cout << endl;
    }
    return 0;
}