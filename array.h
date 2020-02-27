#ifndef __ARRAY_H__
#define __ARRAY_H__

float fa16_dot(const float *a, unsigned int a_stride, const float *b, unsigned int b_stride, unsigned int n);
float fa16_norm(const float *a, unsigned int a_stride, unsigned int n);
void  fa16_unalias(void *dest, void **a, void **b, void *tmp, unsigned size);

#endif /* __ARRAY_H__ */