#ifndef __V3D_H__
#define __V3D_H__

typedef struct {
    float x;
    float y;
    float z;
} v3d_t;

void v3d_add(v3d_t *dest, const v3d_t *a, const v3d_t *b);
void v3d_sub(v3d_t *dest, const v3d_t *a, const v3d_t *b);
void v3d_mul_s(v3d_t *dest, const v3d_t *a, float b);
void v3d_div_s(v3d_t *dest, const v3d_t *a, float b);
// Norm
float v3d_norm(const v3d_t *a);
void v3d_normalize(v3d_t *dest, const v3d_t *a);
// Dot product
float v3d_dot(const v3d_t *a, const v3d_t *b);
// Cross product
void v3d_cross(v3d_t *dest, const v3d_t *a, const v3d_t *b);

#endif /* __V3D_H__ */