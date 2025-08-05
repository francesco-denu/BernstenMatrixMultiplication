#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX 8192*8192 // Possiamo moltiplicare massimo matrici di questa dimensione
typedef unsigned int dtype; // Definisco un tipo in questo modo cosi che sia piu facile cambiarlo dopo

typedef struct {
    int n_rows;        
    int n_cols;         
    dtype entries[MAX]; 
} LOCAL_MATRIX_T;

#define Entry(A, i, j) (*(((A)->entries) + ((A)->n_cols)*(i) + (j))) // Extract elements from LOCAL_MATRIX_T using indices i, j

int read_matrix_from_file(const char *filename, LOCAL_MATRIX_T *matrix) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }

    if (fscanf(file, "%d %d", &matrix->n_rows, &matrix->n_cols) != 2) {
        fclose(file);
        return -1;
    }

    for (int i = 0; i < matrix->n_rows; ++i) {
        for (int j = 0; j < matrix->n_cols; ++j) {
            if (fscanf(file, "%d", &Entry(matrix, i, j)) != 1) {
                fclose(file);
                return -1;
            }
        }
    }

    fclose(file);
    return 0;
}

void multiply_matrices(const LOCAL_MATRIX_T *A, const LOCAL_MATRIX_T *B, LOCAL_MATRIX_T *C) {
    if (A->n_cols != B->n_rows) {
        printf("Error: matrices cannot be multiplied\n");
        return;
    }

    C->n_rows = A->n_rows;
    C->n_cols = B->n_cols;

    for (int i = 0; i < C->n_rows; ++i) {
        for (int j = 0; j < C->n_cols; ++j) {
            Entry(C, i, j) = 0;
            for (int k = 0; k < A->n_cols; ++k) {
                Entry(C, i, j) += Entry(A, i, k) * Entry(B, k, j);
            }
        }
    }
}

void write_matrix_to_file(const char *filename, const LOCAL_MATRIX_T *matrix) {
    FILE *file = fopen(filename, "w+");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    // Write the number of rows and columns
    fprintf(file, "%d %d\n", matrix->n_rows, matrix->n_cols);

    // Write the matrix entries
    for (int i = 0; i < matrix->n_rows; ++i) {
        for (int j = 0; j < matrix->n_cols; ++j) {
            fprintf(file, "%d ", Entry(matrix, i, j));
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

int main(int argc, char *argv[]) {
    clock_t start, end;
    double cpu_time_used;

    if (argc != 4) {
        printf("Usage: %s <matrix1.txt> <matrix2.txt> <output.txt>\n", argv[0]);
        return -1;
    }

    LOCAL_MATRIX_T *A, *B, *C;
    A = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    B = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    C = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));

    // Read matrices from files
    if (read_matrix_from_file(argv[1], A) != 0) {
        printf("Error reading %s\n", argv[1]);
        return -1;
    }

    if (read_matrix_from_file(argv[2], B) != 0) {
        printf("Error reading %s\n", argv[2]);
        return -1;
    }

    start = clock();
    multiply_matrices(A, B, C);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("SEQUENTIAL Algorithm N = %d, Time = %f seconds\n", A->n_rows , cpu_time_used);
    //printf("%f\n", cpu_time_used);
    write_matrix_to_file(argv[3], C);

    return 0;
}
