#include "matrix.h"
#include "array.h"

/****************************
 * Initialization functions *
 ****************************/

void matrix_fill(matrix_t *dest, float value)
{
    int row, column;

    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = value;
        }
    }
}

void matrix_fill_diagonal(matrix_t *dest, float value)
{
    int row;

    matrix_fill(dest, 0);

    for (row = 0; row < dest->rows; row++)
    {
        dest->data[row][row] = value;
    }
}

/*********************************
 * Operations between 2 matrices *
 *********************************/

/* matrix_t *dest must not be matrix_t *a */
void matrix_mul(matrix_t *dest, const matrix_t *a, const matrix_t *b)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp matrix.
    matrix_t tmp;
    array_unalias(dest, (void**)&a, (void**)&b, &tmp, sizeof(matrix_t));
    
    if (a->columns != b->rows)
        return;
    
    dest->rows = a->rows;
    dest->columns = b->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = array_dot(
                &a->data[row][0], 1,
                &b->data[0][column], MATRIX_MAX_SIZE,
                a->columns);
        }
    }

    return;
}

// Multiply transpose of at with b
void matrix_mul_at(matrix_t *dest, const matrix_t *at, const matrix_t *b)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp matrix.
    matrix_t tmp;
    array_unalias(dest, (void**)&at, (void**)&b, &tmp, sizeof(tmp));
    
    if (at->rows != b->rows)
        return;
    
    dest->rows = at->columns;
    dest->columns = b->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = array_dot(
                &at->data[0][row], MATRIX_MAX_SIZE,
                &b->data[0][column], MATRIX_MAX_SIZE,
                at->rows);
        }
    }
}

void matrix_mul_bt(matrix_t *dest, const matrix_t *a, const matrix_t *bt)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp matrix.
    matrix_t tmp;
    array_unalias(dest, (void**)&a, (void**)&bt, &tmp, sizeof(tmp));
    
    if (a->columns != bt->columns)
        return;
    
    dest->rows = a->rows;
    dest->columns = bt->rows;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = array_dot(
                &a->data[row][0], 1,
                &bt->data[column][0], 1,
                a->columns);
        }
    }
}

static inline void matrix_addsub(matrix_t *dest, const matrix_t *a, const matrix_t *b, unsigned int add)
{
    unsigned int row, column;
    float sum;

    if (a->columns != b->columns || a->rows != b->rows)
        return;

    dest->rows = a->rows;
    dest->columns = a->columns;

    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            if (add) {
                sum = a->data[row][column] + b->data[row][column];
            } else {
                sum = a->data[row][column] - b->data[row][column];
            }

            dest->data[row][column] = sum;
        }
    }
}

void matrix_add(matrix_t *dest, const matrix_t *a, const matrix_t *b)
{
    matrix_addsub(dest, a, b, 1);
}

void matrix_sub(matrix_t *dest, const matrix_t *a, const matrix_t *b)
{
    matrix_addsub(dest, a, b, 0);
}

/*********************************
 * Operations on a single matrix *
 *********************************/

void matrix_transpose(matrix_t *dest, const matrix_t *matrix)
{
    int row, column;
    float temp;

    // This code is a bit tricky in order to work
    // in the situation when dest = matrix.
    // Before writing a value in dest, we must copy
    // the corresponding value from matrix to a temporary
    // variable.

    // We actually transpose a n by n square matrix, because
    // that can be done in-place easily. Because matrix_t always
    // allocates a square area even if actual matrix is smaller,
    // this is not a problem.
    unsigned int n = matrix->rows;
    if (matrix->columns > n)
        n = matrix->columns;

    unsigned int rows = matrix->rows;
    dest->rows = matrix->columns;
    dest->columns = rows;

    for (row = 0; row < n; row++)
    {
        for (column = 0; column < row; column++)
        {
            temp = matrix->data[row][column];
            dest->data[row][column] = matrix->data[column][row];
            dest->data[column][row] = temp;
        }

        dest->data[row][row] = matrix->data[row][row];
    }
}

/***************************************
 * Operations of a matrix and a scalar *
 ***************************************/

static inline void matrix_divmul_s(matrix_t *dest, const matrix_t *matrix, float scalar, unsigned int mul)
{
    unsigned int row, column;
    float value;

    dest->rows = matrix->rows;
    dest->columns = matrix->columns;

    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            value = matrix->data[row][column];

            if (mul) {
                value = value * scalar;
            } else {
                value = value / scalar;
            }

            dest->data[row][column] = value;
        }
    }
}

void matrix_mul_s(matrix_t *dest, const matrix_t *matrix, float scalar)
{
    matrix_divmul_s(dest, matrix, scalar, 1);
}

void matrix_div_s(matrix_t *dest, const matrix_t *matrix, float scalar)
{
    matrix_divmul_s(dest, matrix, scalar, 0);
}

/***************************************************
 * Solving linear equations using QR decomposition *
 ***************************************************/

// Takes two columns vectors, v and u, of size n.
// Performs v = v - dot(u, v) * u,
// where dot(u,v) has already been computed
// u is assumed to be an unit vector.
static void subtract_projection(float *v, const float *u, float dot, int n)
{
    float product, diff;
    while (n--)
    {
        product = dot * *u;

        diff = *v - product;

        *v = diff;

        v += MATRIX_MAX_SIZE;
        u += MATRIX_MAX_SIZE;
    }
}

void matrix_qr_decomposition(matrix_t *q, matrix_t *r, const matrix_t *matrix, int reorthogonalize)
{
    int i, j, reorth;
    float dot, norm;

    unsigned int stride = MATRIX_MAX_SIZE;
    unsigned int n = matrix->rows;

    // This uses the modified Gram-Schmidt algorithm.
    // subtract_projection takes advantage of the fact that
    // previous columns have already been normalized.

    // We start with q = matrix
    if (q != matrix)
    {
        *q = *matrix;
    }

    // R is initialized to have square size of cols(A) and zeroed.
    r->columns = matrix->columns;
    r->rows = matrix->columns;
    matrix_fill(r, 0.0f);

    // Now do the actual Gram-Schmidt for the rows.
    for (j = 0; j < q->columns; j++)
    {
        for (reorth = 0; reorth <= reorthogonalize; reorth++)
        {
            for (i = 0; i < j; i++)
            {
                float *v = &q->data[0][j];
                float *u = &q->data[0][i];

                dot = array_dot(v, stride, u, stride, n);
                subtract_projection(v, u, dot, n);

                r->data[i][j] += dot;
            }
        }

        // Normalize the row in q
        norm = array_norm(&q->data[0][j], stride, n);
        r->data[j][j] = norm;

        if (norm < 0.00007629f && norm > -0.00007629f)
        {
            // Nearly zero norm, which means that the row
            // was linearly dependent.
            // Error: singular matrix
            continue;
        }

        for (i = 0; i < n; i++)
        {
            // norm >= v[i] for all i, therefore this division
            // doesn't overflow unless norm approaches 0.
            q->data[i][j] = q->data[i][j] / norm;
        }
    }
}

void matrix_solve(matrix_t *dest, const matrix_t *q, const matrix_t *r, const matrix_t *matrix)
{
    int row, column, variable;

    if (r->columns != r->rows || r->columns != q->columns || r == dest)
    {
        return;
    }

    // Ax=b <=> QRx=b <=> Q'QRx=Q'b <=> Rx=Q'b
    // Q'b is calculated directly and x is then solved row-by-row.
    matrix_mul_at(dest, q, matrix);

    for (column = 0; column < dest->columns; column++)
    {
        for (row = dest->rows - 1; row >= 0; row--)
        {
            float value = dest->data[row][column];

            // Subtract any already solved variables
            for (variable = row + 1; variable < r->columns; variable++)
            {
                float multiplier = r->data[row][variable];
                float known_value = dest->data[variable][column];
                float product = multiplier * known_value;
                value = value - product;
            }

            // Now value = R_ij x_i <=> x_i = value / R_ij
            float divider = r->data[row][row];
            if (divider == 0.0f)
            {
                dest->data[row][column] = 0.0f;
                continue;
            }

            float result = value / divider;
            dest->data[row][column] = result;
        }
    }
}

/**************************
 * Cholesky decomposition *
 **************************/

void matrix_cholesky(matrix_t *dest, const matrix_t *matrix)
{
    // This is the Cholesky–Banachiewicz algorithm.
    // Refer to http://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky.E2.80.93Banachiewicz_and_Cholesky.E2.80.93Crout_algorithms

    int row, column, k;

    if (matrix->rows != matrix->columns)
        return;

    dest->rows = dest->columns = matrix->rows;

    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            if (row == column)
            {
                // Value on the diagonal
                // Ljj = sqrt(Ajj - sum(Ljk^2, k = 1..(j-1))
                float value = matrix->data[row][column];
                for (k = 0; k < column; k++)
                {
                    float Ljk = dest->data[row][k];
                    Ljk = Ljk * Ljk;
                    value = value - Ljk;

                }

                if (value < 0.0f)
                {
                    value = 0.0f;
                }

                dest->data[row][column] = sqrtf(value);
            }
            else if (row < column)
            {
                // Value above diagonal
                dest->data[row][column] = 0.0f;
            }
            else
            {
                // Value below diagonal
                // Lij = 1/Ljj (Aij - sum(Lik Ljk, k = 1..(j-1)))
                float value = matrix->data[row][column];
                for (k = 0; k < column; k++)
                {
                    float Lik = dest->data[row][k];
                    float Ljk = dest->data[column][k];
                    float product = Lik * Ljk;
                    value = value - product;
                }
                float Ljj = dest->data[column][column];
                value = value / Ljj;
                dest->data[row][column] = value;
            }
        }
    }
}

/***********************************
 * Lower-triangular matrix inverse *
 **********************************/

void matrix_invert_lt(matrix_t *dest, const matrix_t *matrix)
{
    // This is port of the algorithm as found in the Efficient Java Matrix Library
    // https://code.google.com/p/efficient-java-matrix-library

    int i, j, k;
    const unsigned int n = matrix->rows;

    // If dest and input matrices alias, we have to use a temp matrix.
    matrix_t tmp;
    array_unalias(dest, (void **)&matrix, (void **)&matrix, &tmp, sizeof(tmp));

    // TODO reorder these operations to avoid cache misses

    // inverts the lower triangular system and saves the result
    // in the upper triangle to minimize cache misses
    for (i = 0; i < n; ++i)
    {
        const float el_ii = matrix->data[i][i];
        for (j = 0; j <= i; ++j)
        {
            float sum = (i == j) ? 1.0f : 0.0f;
            for (k = i - 1; k >= j; --k)
            {
                sum = sum - (matrix->data[i][k] * dest->data[j][k]);
            }
            dest->data[j][i] = sum / el_ii;
        }
    }
    // solve the system and handle the previous solution being in the upper triangle
    // takes advantage of symmetry
    for (i = n - 1; i >= 0; --i)
    {
        const float el_ii = matrix->data[i][i];
        for (j = 0; j <= i; ++j)
        {
            float sum = (i < j) ? 0.0f : dest->data[j][i];
            for (k = i + 1; k < n; ++k)
            {
                sum = sum - (matrix->data[k][i] * dest->data[j][k]);
            }
            dest->data[i][j] = dest->data[j][i] = sum / el_ii;
        }
    }
}
