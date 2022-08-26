#include <omp.h>
#include <cstdio>
using namespace std;

int main() {
#pragma omp parallel
    { printf("Hello, World\n"); }
    printf("Hello\n");
    return 0;
}
