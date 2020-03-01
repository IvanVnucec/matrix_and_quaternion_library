#include "quat.h"
#include "v3d.h"

#include <math.h>

// Conjugate of quaternion
void quat_conj(quat_t *dest, const quat_t *q)
{
    dest->a =   q->a;
    dest->b = - q->b;
    dest->c = - q->c;
    dest->d = - q->d;
}

// Multiply two quaternions, dest = a * b.
void quat_mul(quat_t *dest, const quat_t *q, const quat_t *r) {

    dest->a = (q->a * r->a) - (q->b * r->b) - (q->c * r->c) - (q->d * r->d);
    dest->b = (q->a * r->b) + (q->b * r->a) + (q->c * r->d) - (q->d * r->c);
    dest->c = (q->a * r->c) - (q->b * r->d) + (q->c * r->a) + (q->d * r->b);
    dest->d = (q->a * r->d) + (q->b * r->c) - (q->c * r->b) + (q->d * r->a);

    return;
}

void quat_add(quat_t *dest, const quat_t *q, const quat_t *r)
{
    dest->a = q->a + r->a;
    dest->b = q->b + r->b;
    dest->c = q->c + r->c;
    dest->d = q->d + r->d;
}

// Multiply quaternion by scalar
void quat_mul_s(quat_t *dest, const quat_t *q, float s)
{
    dest->a = (q->a * s);
    dest->b = (q->b * s);
    dest->c = (q->c * s);
    dest->d = (q->d * s);
}

// Divide quaternion by scalar
void quat_div_s(quat_t *dest, const quat_t *q, float s)
{
    dest->a = (q->a / s);
    dest->b = (q->b / s);
    dest->c = (q->c / s);
    dest->d = (q->d / s);
}

float quat_dot(const quat_t *q, const quat_t *r)
{
    return q->a * r->a + q->b * r->b + q->c * r->c + q->d * r->d;    
}

// Quaternion norm
float quat_norm(const quat_t *q)
{
    return sqrtf(quat_dot(q, q));
}

// Normalize quaternion
void quat_normalize(quat_t *dest, const quat_t *q)
{
    quat_div_s(dest, q, quat_norm(q));
}

// Quaternion power
// See https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
void quat_pow(quat_t *dest, const quat_t *q, float power)
{   
    const float norm = quat_norm(q);
    const float norm_pow = powf(norm, power);
    const float angle = acosf(q->a / norm);

    dest->a = norm_pow * cosf(power * angle);
    dest->b = q->b / norm / sinf(angle) * norm_pow * sinf(power * angle);
    dest->c = q->c / norm / sinf(angle) * norm_pow * sinf(power * angle);
    dest->d = q->d / norm / sinf(angle) * norm_pow * sinf(power * angle);

    return;
}

// Weighted average
// See http://www.acsu.buffalo.edu/~johnc/ave_sfm07.pdf
void quat_avg(quat_t *dest, const quat_t *q1, const quat_t *q2, float weight)
{
    // z = sqrt((w1 - w2)^2 + 4 w1 w2 (q1' q2)^2
    // <=>
    // z = sqrt((2 w1 - 1)^2 + 4 w1 (1 - w1) (q1' q2)^2)
    const float dot = quat_dot(q1, q2);
    const float z = sqrtf( (2 * weight - 1.0f) * (2.0f * weight - 1.0f) 
                    + 4.0f * weight * (1.0f - weight) * dot * dot );
    
    // q = 2 * w1 * (q1' q2) q1 + (w2 - w1 + z) q2
    // <=>
    // q = 2 * w1 * (q1' q2) q1 + (1 - 2 * w1 + z) q2
    quat_t tmp1;
    quat_mul_s(&tmp1, q1, 2.0f * weight * dot);
    
    quat_t tmp2;
    quat_mul_s(&tmp2, q2, 1.0f - 2.0f * weight + z);
    
    quat_add(dest, &tmp1, &tmp2);
    quat_normalize(dest, dest);
}


void quat_from_axis_angle(quat_t *dest, const v3d_t *axis, float angle)
{
    const float scale = sinf(angle / 2.0f);
    
    dest->a = cosf(angle);
    dest->b = axis->x * scale;
    dest->c = axis->y * scale;
    dest->d = axis->z * scale;
}


// Unit quaternion to rotation matrix
void quat_to_matrix(matrix_t *dest, const quat_t *q)
{
    dest->rows = dest->columns = 3;

    dest->data[0][0] = 1.0f - 2.0f * (q->c)*(q->c) + (q->d)*(q->d);
    dest->data[1][1] = 1.0f - 2.0f * (q->b)*(q->b) + (q->d)*(q->d);
    dest->data[2][2] = 1.0f - 2.0f * (q->b)*(q->b) + (q->c)*(q->c);
    
    dest->data[1][0] = 2.0f * ((q->b)*(q->c) + (q->a)*(q->d));
    dest->data[0][1] = 2.0f * ((q->b)*(q->c) - (q->a)*(q->d));

    dest->data[2][0] = 2.0f * ((q->b)*(q->d) - (q->a)*(q->c));
    dest->data[0][2] = 2.0f * ((q->b)*(q->d) + (q->a)*(q->c));
    
    dest->data[2][1] = 2.0f * ((q->c)*(q->d) + (q->a)*(q->b));
    dest->data[1][2] = 2.0f * ((q->c)*(q->d) - (q->a)*(q->b));
}


void quat_rotate(v3d_t *dest, const quat_t *q, const v3d_t *v)
{
    quat_t vector, q_conj;
    
    quat_from_v3d(&vector, v, 0);
    quat_conj(&q_conj, q);
    
    quat_mul(&vector, q, &vector);
    quat_mul(&vector, &vector, &q_conj);
    
    quat_to_v3d(dest, &vector);
}
