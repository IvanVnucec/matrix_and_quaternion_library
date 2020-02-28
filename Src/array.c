#include "array.h"

#include <math.h>
#include <string.h>

float fa16_dot(const float *a, unsigned int a_stride, const float *b, unsigned int b_stride, unsigned int n) {
    float sum = 0.0;

    while(n--) {
        sum += *a * *b;

        a += a_stride;
        b += b_stride;
    }

    return sum;
}


float fa16_norm(const float *a, unsigned int a_stride, unsigned int n)
{
    float sum = 0.0f;
    
    while (n--)
    {
        if (*a != 0)
        {
            sum += (*a) * (*a);
        }
        
        a += a_stride;
    }
    
    return sqrtf(sum);
}


void fa16_unalias(void *dest, void **a, void **b, void *tmp, unsigned size)
{
    if (dest == *a)
    {
        memcpy(tmp, *a, size);
        *a = tmp;
        
        if (dest == *b)
            *b = tmp;
    }
    else if (dest == *b)
    {
        memcpy(tmp, *b, size);
        *b = tmp;
    }
}