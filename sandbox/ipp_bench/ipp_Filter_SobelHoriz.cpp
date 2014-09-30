#include <stdio.h>

#include <ippcore.h>
#include <ippvm.h>
#include <ipps.h>
#include <ippi.h>

#include "timing.h"

# define nChannels 3

int main ()
{
    timeval t1, t2;

    Ipp8u pSrc[8*8]={
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 3, 3, 3, 3, 3, 3, 4,
    4, 3, 2, 2, 2, 2, 3, 4,
    4, 3, 2, 1, 1, 2, 3, 4,
    4, 3, 2, 1, 1, 2, 3, 4,
    4, 3, 2, 2, 2, 2, 3, 4,
    4, 3, 3, 3, 3, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4
    };

    Ipp16s pDst[8*8];
    IppiSize Roi = { 8, 8 };

    timer__(&t1); // start clock

    // apply second order sobel vertical filter
    ippiFilterSobelVertSecond_8u16s_C1R ( pSrc, 8, pDst, 8*sizeof( Ipp16s ),
        Roi, ippMskSize3x3 );

    __timer(&t1, &t2, "conversion"); // end clock

    // print dst array
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            printf("%d, ", pDst[i*8 + j]);
        }
        printf("\n");
    }

    return 0; 
}
