#ifndef __QUAT_H__
#define __QUAT_H__

#include "v3d.h"
#include "matrix.h"

typedef struct {
    float a;
    float b;
    float c;
    float d;
} quat_t;

void  quat_conj (quat_t *dest, const quat_t *q);
void  quat_mul  (quat_t *dest, const quat_t *q, const quat_t *r);
void  quat_add  (quat_t *dest, const quat_t *q, const quat_t *r);
void  quat_mul_s(quat_t *dest, const quat_t *q, float s);
void  quat_div_s(quat_t *dest, const quat_t *q, float s);
float quat_dot  (const quat_t *q, const quat_t *r);
float quat_norm (const quat_t *q);
void  quat_normalize(quat_t *dest, const quat_t *q);
void  quat_pow  (quat_t *dest, const quat_t *q, float power);
void  quat_avg  (quat_t *dest, const quat_t *q1, const quat_t *q2, float weight);
void  quat_from_axis_angle(quat_t *dest, const v3d_t *axis, float angle);
void  quat_to_matrix(matrix_t *dest, const quat_t *q);
void  quat_rotate   (v3d_t *dest, const quat_t *q, const v3d_t *v);

static inline void quat_from_v3d(quat_t *q, const v3d_t *v, float a)
{
    q->a = a;
    q->b = v->x;
    q->c = v->y;
    q->d = v->z;
}

static inline void quat_to_v3d(v3d_t *v, const quat_t *q)
{
    v->x = q->b;
    v->y = q->c;
    v->z = q->d;
}

#endif /* __QUAT_H__ */