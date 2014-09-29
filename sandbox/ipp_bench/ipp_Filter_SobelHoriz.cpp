#include <stdio.h>

#include <ippcore.h>
#include <ippvm.h>
#include <ipps.h>
#include <ippi.h>

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

# define nChannels 3

int main ()
{
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

    ippiFilterSobelVertSecond_8u16s_C1R ( pSrc, 8, pDst, 8*sizeof( Ipp16s ),
        Roi, ippMskSize3x3 );

    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            printf("%d, ", pDst[i*8 + j]);
        }
        printf("\n");
    }

    return 0; 
}
