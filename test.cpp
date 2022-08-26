#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <atomic>
#define THREADS 4
#define N 16
#define CHUNK 100
using namespace std;

atomic_int total(0);

int main() {
    clock_t start, end;
    start = clock();

    int i;
#pragma omp parallel for schedule(guided) num_threads(THREADS)
    for (i = 0; i < N; i++) {
        /* wait for i seconds */
        sleep(i);
        total++;
        printf("Thread %d has completed iteration %d.\n", omp_get_thread_num(), i);
    }
    /* all threads done */
    printf("All done!\n");
    printf("Total %u\n", total.load());
    end = clock();
    printf("Run Time : %.10lf\n", (double)(end - start) / 1000.0);
    return 0;
}