#include <omp.h>
#include <cstdio>
using namespace std;

// g++ --std=c++20 -Xpreprocessor -fopenmp -I/opt/homebrew/Cellar/libomp/14.0.6/include -L/opt/homebrew/Cellar/llvm/14.0.6_1/lib -lomp test.cpp 


int main() {
#pragma omp parallel
    { printf("Hello, World\n"); }
    printf("Hello\n");
    return 0;
}
