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
    bool* a = new bool[10];
    for (int i = 0; i < 10; i++) {
        a[i] = true;
    }
    bool* b = new bool[10];
    for (int i = 0; i < 10; i++) {
        b[i] = a[i];
    }

    cout << " a : " << a << endl;
    cout << " b : " << b << endl;

    a[0] = false; // 
    swap(a, b);
    cout << " a : " << a << endl;
    cout << " b : " << b << endl;

    for (int i = 0; i < 10; i++) {
        cout << a[i] << " ";
    }
    cout << endl;
    for (int i = 0; i < 10; i++) {
        cout << b[i] << " ";
    }
    cout << endl;

    
    return 0;
}