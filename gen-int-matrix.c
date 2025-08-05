#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <filename> <n> <m>\n", argv[0]);
        return 1;
    }

    char *filename = argv[1];
    int n = atoi(argv[2]);
    int m = atoi(argv[3]);

    if (n <= 0 || m <= 0) {
        fprintf(stderr, "Both n and m must be positive integers.\n");
        return 1;
    }

    FILE *file = fopen(filename, "w+");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    fprintf(file, "%d %d\n", n, m);
    
    int count = 1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            //fprintf(file, "%d ", count++);
            int random_value = (rand() % 2) + 1;
            fprintf(file, "%d ", random_value);
        }
        fprintf(file, "\n");
    }

    fclose(file);

    printf("Matrix written to %s\n", filename);
    return 0;
}
