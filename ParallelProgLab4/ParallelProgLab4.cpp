//#define _CRT_SECURE_NO_WARNINGS
//
//#include <stdio.h>
//#include <corecrt_math.h>
//#include <stdlib.h>
//#include <ctime>
//#include <windows.h>
//#include "Header.h"
//#include "../../../../../../Program Files (x86)/Microsoft SDKs/MPI/Include/mpi.h"
//#include <locale.h>
//#include <corecrt_math_defines.h>
//#pragma comment(lib, "../../../../../../Program Files (x86)/Microsoft SDKs/MPI/Lib/x86/msmpi.lib")
//
//double f(double x) {
//    return sin(x); // Подынтегральная функция
//}
//
//double integrate(double a, double b, int n) {
//    double h = (b - a) / n;
//    double sum = 0.0;
//
//    for (int i = 0; i < n; i++) {
//        double x = a + (i + 0.5) * h; // Средняя точка
//        sum += f(x);
//    }
//
//    return sum * h;
//}
//
//int main(int argc, char** argv) {
//    #ifdef _WIN32
//        SetConsoleCP(CP_UTF8);
//        SetConsoleOutputCP(CP_UTF8);
//    #endif
//    int rank, size;
//    double a = 0.0; // Левая граница
//    double b = M_PI; // Правая граница (π)
//    int total_intervals = 1000000; // Общее количество интервалов
//    int intervals_per_process;
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    // Расчет количества интервалов для каждого процесса
//    intervals_per_process = total_intervals / size;
//
//    // Определение локальных границ интегрирования
//    double local_a = a + rank * intervals_per_process * (b - a) / total_intervals;
//    double local_b = local_a + intervals_per_process * (b - a) / total_intervals;
//
//    // Начало замера времени
//    double start_time = MPI_Wtime();
//
//    printf("1");
//
//    // Локальное интегрирование
//    double local_result = integrate(local_a, local_b, intervals_per_process);
//
//    // Сбор результатов на нулевом процессе
//    double* results = NULL;
//
//    if (rank == 0) {
//        results = (double*)malloc(size * sizeof(double)); // Выделение памяти для результатов
//    }
//
//    MPI_Gather(&local_result, 1, MPI_DOUBLE, results, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    // Конец замера времени
//    double end_time = MPI_Wtime();
//
//    if (rank == 0) {
//        // Суммирование результатов на нулевом процессе
//        double global_result = 0.0;
//        for (int i = 0; i < size; i++) {
//            global_result += results[i];
//        }
//
//        printf("Result: %f\n", global_result);
//        printf("Time: %f seconds\n", end_time - start_time);
//
//        free(results); // Освобождение памяти
//    }
//
//    MPI_Finalize();
//
//    return 0;
//}
