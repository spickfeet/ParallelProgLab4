#define _CRT_SECURE_NO_WARNINGS
#include <corecrt_math.h>
#include <stdlib.h>
#include <ctime>
#include <windows.h>
#include "Header.h"
#include "../../../../../../Program Files (x86)/Microsoft SDKs/MPI/Include/mpi.h"
#include <locale.h>
#include <corecrt_math_defines.h>
#include <stdio.h>
#pragma comment(lib, "../../../../../../Program Files (x86)/Microsoft SDKs/MPI/Lib/x86/msmpi.lib")

//const int ROWS = 4; // Количество строк матрицы
//const int COLS = 4;
//
//void generate_matrix_and_vector(double matrix[ROWS][COLS], double vector[], int fill_with_ones) {
//if (fill_with_ones) {
//        for (int i = 0; i < ROWS; i++)
//            for (int j = 0; j < COLS; j++)
//                matrix[i][j] = 1.0;
//
//        for (int j = 0; j < COLS; j++)
//            vector[j] = 1.0;
//    }
//    else {
//        srand(time(NULL));
//        for (int i = 0; i < ROWS; i++)
//            for (int j = 0; j < COLS; j++)
//                matrix[i][j] = rand() % 10; // случайные числа от 0 до 9
//
//        for (int j = 0; j < COLS; j++)
//            vector[j] = rand() % 10;
//    }
//}
//
//void print_vector(int size, double vector[ROWS]) {
//    printf("Resulting vector:\n");
//    for (int i = 0; i < size; i++) {
//        printf("%f ", vector[i]);
//    }
//    printf("\n");
//}
//
//int main(int argc, char** argv) {
//    int rank, size;
//    const int ROWS = 4; // Количество строк матрицы
//    const int COLS = 4; // Количество столбцов матрицы
//    double matrix[ROWS][COLS]; // Исходная матрица
//    double vector[COLS]; // Исходный вектор
//    double local_result[ROWS]; // Локальный результат для каждого процесса
//    double global_result[ROWS]; // Глобальный результат для нулевого процесса
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    if (rank == 0) {
//        generate_matrix_and_vector(matrix, vector, 1); // Заполнение единицами
//        printf("Matrix:\n");
//        for (int i = 0; i < ROWS; i++) {
//            for (int j = 0; j < COLS; j++) {
//                printf("%f ", matrix[i][j]);
//            }
//            printf("\n");
//        }
//        printf("Vector:\n");
//        for (int j = 0; j < COLS; j++) {
//            printf("%f ", vector[j]);
//        }
//        printf("\n");
//    }
//
//    // Рассылка строк матрицы другим процессам
//    MPI_Bcast(vector, COLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    // Определение количества строк для каждого процесса
//    int rows_per_process = ROWS / size;
//
//    // Создание производного типа данных для строки матрицы
//    MPI_Datatype row_type;
//    MPI_Type_create_subarray(2, { ROWS, COLS }, (int[]) { 1, COLS }, (int[]) { rank* rows_per_process, 0 }, MPI_ORDER_C, MPI_DOUBLE, & row_type);
//    MPI_Type_commit(&row_type);
//
//    // Локальная часть матрицы для текущего процесса
//    double local_matrix[rows_per_process][COLS];
//
//    // Рассылка строк матрицы по процессам
//    MPI_Scatter(matrix, rows_per_process * COLS, MPI_DOUBLE,
//        local_matrix, rows_per_process * COLS, MPI_DOUBLE,
//        0, MPI_COMM_WORLD);
//
//    // Умножение локальной части матрицы на вектор
//    for (int i = 0; i < rows_per_process; i++) {
//        local_result[i] = 0.0;
//        for (int j = 0; j < COLS; j++) {
//            local_result[i] += local_matrix[i][j] * vector[j];
//        }
//    }
//
//    // Сбор результатов на нулевом процессе
//    MPI_Gather(local_result, rows_per_process, MPI_DOUBLE,
//        global_result, rows_per_process, MPI_DOUBLE,
//        0, MPI_COMM_WORLD);
//
//    if (rank == 0) {
//        print_vector(ROWS, global_result); // Печать результата
//    }
//
//    // Освобождение ресурсов и завершение работы с MPI
//    MPI_Type_free(&row_type);
//    MPI_Finalize();
//
//    return 0;
//}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void multiply_block(int* block, int* vector, int* result, int block_size, int N) {
    for (int i = 0; i < block_size; i++) {
        result[i] = 0;
        for (int j = 0; j < N; j++) {
            result[i] += block[i * N + j] * vector[j];
        }
    }
}

int main(int argc, char* argv[]) {
    int rank, size;
    const int M = 4, N = 4;  // Размеры матрицы и вектора
    int* matrix = NULL;
    int* vector = NULL;
    int* result = NULL;
    MPI_Datatype block_type, resized_block_type;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Размер блока для каждого процесса
    int block_size = M / size;
    if (M % size != 0) {
        if (rank == 0) {
            printf("Ошибка: матрица не делится на равные блоки по строкам!\n");
        }
        MPI_Finalize();
        return -1;
    }

    if (rank == 0) {
        matrix = (int*)malloc(M * N * sizeof(int));
        vector = (int*)malloc(N * sizeof(int));
        result = (int*)malloc(M * sizeof(int));

        // Генерация данных
        srand(time(NULL));
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                matrix[i * N + j] = i * N + j;
            }
        }
        for (int i = 0; i < N; i++) {
            vector[i] = i;
        }
    }

    double start_time = MPI_Wtime();

    // Рассылка вектора всем процессам
    if (rank != 0) {
        vector = (int*)malloc(N * sizeof(int));
    }
    MPI_Bcast(vector, N, MPI_INT, 0, MPI_COMM_WORLD);

    // Создание пользовательского типа данных
    MPI_Type_vector(block_size, N, N, MPI_INT, &block_type);
    MPI_Type_create_resized(block_type, 0, block_size * N * sizeof(int), &resized_block_type);
    MPI_Type_commit(&resized_block_type);

    // Распределение блоков матрицы
    int* recv_block = (int*)malloc(block_size * N * sizeof(int));
    MPI_Scatter(matrix, 1, resized_block_type, recv_block, 1, resized_block_type, 0, MPI_COMM_WORLD);

    // Выделение памяти для локального результата
    int* local_result = (int*)malloc(block_size * sizeof(int));

    // Умножение блока на вектор
    multiply_block(recv_block, vector, local_result, block_size, N);

    // Сбор результатов
    MPI_Gather(local_result, block_size, MPI_INT, result, block_size, MPI_INT, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();
    // Печать результатов
    if (rank == 0) {
        printf("Matrix:\n");
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                printf("%d ", matrix[i * N + j]);
            }
            printf("\n");
        }
        printf("Vector:\n");
        for (int i = 0; i < M; i++) {
            printf("%d\n", vector[i]);
        }
        printf("Result:\n");
        for (int i = 0; i < M; i++) {
            printf("%d\n", result[i]);
        }

        printf("Time: %f seconds\n", end_time - start_time);
        free(matrix);
        free(vector);
        free(result);
    }

    free(recv_block);
    free(local_result);

    MPI_Type_free(&block_type);
    MPI_Type_free(&resized_block_type);

    MPI_Finalize();
    return 0;
}