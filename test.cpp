#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <atomic>
#define THREADS 16
#define N 10000
#define CHUNK 100
using namespace std;

atomic_int total(0);

int main() {
    clock_t start, end;
    start = clock();

    int i;
#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
    for (i = 0; i < N; i++) {
#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
        for (int j = 0; j < N; j++) {
            total++;
        }
        /* wait for i seconds */
    }
    /* all threads done */
    printf("All done!\n");
    printf("Total %u\n", total.load());
    end = clock();
    printf("Run Time : %.10lf\n", (double)(end - start) / CLOCKS_PER_SEC);
    return 0;
}