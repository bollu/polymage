#include <stdio.h>

void gray(int C, int R, const void * img, void *result)
{
    unsigned char *imgc = (unsigned char *)img;
    unsigned char *resultc = (unsigned char *)result;

    int i, j;

    //printf("%d %d %d %d\n", C, R, img, result);
    #pragma omp parallel for
    for(i = 0; i < R; i++){
        for(j = 0; j < C; j++){
            resultc[i*C + j] =
                               imgc[i*3*C + j*3 + 0] * 0.299f +
                               imgc[i*3*C + j*3 + 1] * 0.587f +
                               imgc[i*3*C + j*3 + 2] * 0.114f;
        }
    }
}
