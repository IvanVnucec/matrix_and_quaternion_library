#ifndef __MATRIX_H__
#define __MATRIX_H__

#define MATRIX_MAX_SIZE 8

typedef struct {
    unsigned int rows;
    unsigned int columns;

    float data[MATRIX_MAX_SIZE][MATRIX_MAX_SIZE];
} matrix_t;

#include "matrix.h"

#include <string.h>
#include <math.h>


/****************************
 * Initialization functions *
 ****************************/

void matrix_fill(matrix_t *dest, float value);

void matrix_fill_diagonal(matrix_t *dest, float value);

/*********************************
 * Operations between 2 matrices *
 *********************************/

/* matrix_t *dest must not be matrix_t *a */
void matrix_mul(matrix_t *dest, const matrix_t *a, const matrix_t *b);

// Multiply transpose of at with b
void matrix_mul_at(matrix_t *dest, const matrix_t *at, const matrix_t *b);

void matrix_mul_bt(matrix_t *dest, const matrix_t *a, const matrix_t *bt);

void matrix_add(matrix_t *dest, const matrix_t *a, const matrix_t *b);

void matrix_sub(matrix_t *dest, const matrix_t *a, const matrix_t *b);

/*********************************
 * Operations on a single matrix *
 *********************************/

void matrix_transpose(matrix_t *dest, const matrix_t *matrix);

/***************************************
 * Operations of a matrix and a scalar *
 ***************************************/

void matrix_mul_s(matrix_t *dest, const matrix_t *matrix, float scalar);

void matrix_div_s(matrix_t *dest, const matrix_t *matrix, float scalar);

/***************************************************
 * Solving linear equations using QR decomposition *
 ***************************************************/

void matrix_qr_decomposition(matrix_t *q, matrix_t *r, const matrix_t *matrix, int reorthogonalize);

void matrix_solve(matrix_t *dest, const matrix_t *q, const matrix_t *r, const matrix_t *matrix);
/**************************
 * Cholesky decomposition *
 **************************/

void matrix_cholesky(matrix_t *dest, const matrix_t *matrix);

/***********************************
 * Lower-triangular matrix inverse *
 **********************************/

void matrix_invert_lt(matrix_t *dest, const matrix_t *matrix);

#endif /* __MATRIX_H__ */