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

void matrix_vector_multiply(double* matrix, double* vector, double* result, int rows, int cols, int rank, int size) {
    // Каждый процесс вычисляет часть результата
    int rows_per_process = rows / size;
    int start_row = rank * rows_per_process;
    int end_row = (rank == size - 1) ? rows : start_row + rows_per_process;

    // Локальные вычисления
    for (int i = start_row; i < end_row; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i * cols + j] * vector[j];
        }
    }
}

int main(int argc, char** argv) {
    srand(time(0));
    int rank, size;
    int rows = 3;  // Примерный размер матрицы
    int cols = 3;
    double* matrix = NULL, * vector = NULL, * result = NULL;

    // Инициализация MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Генерация данных
    if (rank == 0) {
        matrix = (double*)malloc(rows * cols * sizeof(double));
        vector = (double*)malloc(cols * sizeof(double));
        result = (double*)malloc(rows * sizeof(double));

        // Заполнение матрицы и вектора единицами
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i * cols + j] = 1.0 * (rand() % 10);  // Заполнение единицами
            }
        }
        for (int j = 0; j < cols; j++) {
            vector[j] = 1.0 * (rand() % 10);  // Заполнение единицами
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                printf("%f ", matrix[i * cols + j]);  // Заполнение единицами
            }
            printf("\n");
        }
        printf("\n");
        for (int j = 0; j < cols; j++) {
            printf("%f ", vector[j]);  // Заполнение единицами
        }
        printf("\n");
    }

    // Рассылка данных
    double* local_matrix = (double*)malloc(rows / size * cols * sizeof(double));
    double* local_result = (double*)malloc(rows / size * sizeof(double));

    // Разделение матрицы по строкам
    MPI_Scatter(matrix, rows / size * cols, MPI_DOUBLE, local_matrix, rows / size * cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Разделение вектора по всем процессам (вектор один и тот же для всех)
    MPI_Bcast(vector, cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Умножение подматрицы на вектор
    matrix_vector_multiply(local_matrix, vector, local_result, rows, cols, rank, size);

    // Сбор результата в процессе 0
    MPI_Gather(local_result, rows / size, MPI_DOUBLE, result, rows / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Вывод результата
    if (rank == 0) {
        for (int i = 0; i < rows; i++) {
            printf("result[%d] = %f\n", i, result[i]);
        }
        free(matrix);
        free(vector);
        free(result);
    }

    // Очистка
    free(local_matrix);
    free(local_result);

    MPI_Finalize();
    return 0;
}