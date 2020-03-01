#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "v3d.h"
#include "quat.h"
#include "matrix.h"

void print_quaternion(quat_t q) {
    printf("q = (%f %f %f %f)\n", q.a, q.b, q.c, q.d);
}

void print_matrix(matrix_t m) {
    for (int i=0; i<m.rows; i++) {
        for (int j=0; j<m.columns; j++) {
            printf("%f ", m.data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(v3d_t v) {
    printf("v = (%f %f %f)\n", v.x, v.y, v.z);
}

int main(void) {
    //v3d_t vector = {1.0f, 2.0f, 3.0f};
    //quat_t quaternion = {1.0f, 2.0f, 3.0f, 4.0f};
    matrix_t matrix1 = {3, 3, {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}};
    matrix_t matrix2 = {3, 3, {{1, 0, -1}, {0, 1, 0}, {1, 2, 3}}};
    matrix_t matrix3;

    print_matrix(matrix1);
    print_matrix(matrix2);

    matrix_horzcat(&matrix3, &matrix1, &matrix2);
    print_matrix(matrix3);

    matrix_vertcat(&matrix3, &matrix1, &matrix2);
    print_matrix(matrix3);
    
    return 0;
}