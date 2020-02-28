#include "v3d.h"
#include "array.h"

void v3d_add(v3d_t *dest, const v3d_t *a, const v3d_t *b)
{
        dest->x = a->x + b->x;
        dest->y = a->y + b->y;
        dest->z = a->z + b->z;
}

void v3d_sub(v3d_t *dest, const v3d_t *a, const v3d_t *b)
{
        dest->x = a->x - b->x;
        dest->y = a->y - b->y;
        dest->z = a->z - b->z;
}

void v3d_mul_s(v3d_t *dest, const v3d_t *a, float b)
{
        dest->x = a->x * b;
        dest->y = a->y * b;
        dest->z = a->z * b;
}

void v3d_div_s(v3d_t *dest, const v3d_t *a, float b)
{
        dest->x = a->x / b;
        dest->y = a->y / b;
        dest->z = a->z / b;
}

// Norm
float v3d_norm(const v3d_t *a)
{
    return fa16_norm(&a->x, &a->y - &a->x, 3);
}

void v3d_normalize(v3d_t *dest, const v3d_t *a)
{
    v3d_div_s(dest, a, v3d_norm(a));
}

// Dot product
float v3d_dot(const v3d_t *a, const v3d_t *b)
{
    return fa16_dot(&a->x, &a->y - &a->x, &b->x, &b->y - &b->x, 3);
}

// Cross product
void v3d_cross(v3d_t *dest, const v3d_t *a, const v3d_t *b)
{
    v3d_t tmp;
    fa16_unalias(dest, (void**)&a, (void**)&b, &tmp, sizeof(tmp));
    
    dest->x = a->y*b->z - a->z*b->y;
    dest->y = a->z*b->x - a->x*b->z;
    dest->z = a->x*b->y - a->y*b->x;
}