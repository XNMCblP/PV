#include <iostream>
#include "cstdlib"
#include <omp.h>

int main() {
    int Q = 14000000;
    int Par = 1;
    double t1, t2;
    srand(222);
    int N = 5;
    double* a = new double[N * N];
    double* b = new double[N * N];
    double* c = new double[N * N];
    for (int i = 0; i < N * N; i++) {
        a[i] = rand() * 1.0 / RAND_MAX - 0.5;
        b[i] = rand() * 1.0 / RAND_MAX - 0.5;
        c[i] = 0;
    }
    omp_set_num_threads(4);
    t1 = omp_get_wtime();
    if (Par == 1) {
#pragma omp parallel shared (c)
        {
        for (int q = 0; q < Q; q++) {
            double* mas = new double[N * N];
            for (int i = 0; i < N * N; i++) {
                mas[i] = 0;
            }
#pragma omp for
            for (int k = 0; k < N; k++) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        mas[i * N + j] = a[i * N + k] * b[k * N + i];
                    }
                }
            }
#pragma omp critical
                for (int i = 0; i < N * N; i++) {
                    c[i] += mas[i];
                }
            }
        }
    }
    else if (Par == 2) {
#pragma omp parallel shared (c)
        {
            for (int q = 0; q < Q; q++) {
                double* mas = new double[N * N];
                for (int i = 0; i < N * N; i++) {
                    mas[i] = 0;
                }
#pragma omp for
                for (int k = 0; k < N; k++) {
                    for (int j = 0; j < N; j++) {
                        for (int i = 0; i < N; i++) {
                            mas[i * N + j] = a[i * N + k] * b[k * N + i];
                        }
                    }
                }
#pragma omp critical
                for (int i = 0; i < N * N; i++) {
                    c[i] += mas[i];
                }
            }
        }
    }
    else if (Par == 3) {
        for (int q = 0; q < Q; q++) {
            for(int k = 0; k < N; k++){
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        c[i * N + j] = a[i * N + k] * b[k * N + i];
                    }
                }
            }
        }
    }
    else {
#pragma omp parallel shared(c)
        {
            for (int q = 0; q < Q; q++) {
                double* mas = new double[N * N];
                for (int i = 0; i < N * N; i++) {
                    mas[i] = 0;
                }
#pragma omp for
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        for (int k = 0; k < N; k++) {
                            mas[i * N + j] = a[i * N + k] * b[k * N + i];
                        }
                    }
                }
#pragma omp critical
                for (int i = 0; i < N * N; i++) {
                    c[i] += mas[i];
                }
            }
        }
    }
    t2 = omp_get_wtime();
    std::cout << "Time:" << t2 - t1 << std::endl;
}