#include <stdio.h>
#include <stdlib.h>

int* convolve(int* u, int u_size, int* v, int v_size) {
    // Set Results size and elems to zero
    int result_size = u_size + v_size - 1;
    int* result = (int*) calloc(result_size, sizeof(int));
    
    for (int i = 0; i < result_size; i++) {
        result[i] = 0;
        for (int j = 0; j < u_size; j++) {
            if (i - j >= 0 && i - j < v_size) {
                result[i] += u[j] * v[i - j];
            }
        }
    }
    return result;
}

int main() {
    int u[] = {1, 2, 3, 4, 5};
    int v[] = {1, 1, 0};
    int u_size = sizeof(u) / sizeof(u[0]);
    int v_size = sizeof(v) / sizeof(v[0]);

    int* result = convolve(u, u_size, v, v_size);
    int result_size = u_size + v_size - 1;

    printf("Convolution result: ");
    for (int i = 0; i < result_size; i++) {
        printf("%d ", result[i]);
    }
    printf("\n");

    free(result);
    return 0;
}
