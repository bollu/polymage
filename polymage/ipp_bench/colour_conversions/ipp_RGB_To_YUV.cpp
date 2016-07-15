#include <stdio.h>
#include <sys/time.h>

#include <ippcore.h>
#include <ippvm.h>
#include <ipps.h>
#include <ippi.h>
#include <ippcc.h>

#include "timing.h"

# define nChannels 3

int main ()
{
    timeval t1, t2;

    Ipp8u src [3*3*nChannels] = {
        255, 0, 0, 255, 0, 0, 255, 0, 0,
        0, 255, 0, 0, 255, 0, 0, 255, 0,
        0, 0, 255, 0, 0, 255, 0, 0, 255};

    Ipp8u dst [3*3* nChannels ];

    IppiSize roiSize = { 3, 3 };

    int srcStep = 3 * nChannels ;
    int dstStep = 3 * nChannels ; 

    timer__(&t1); // start clock
    // conversion from RGB to YUV
    ippiRGBToYUV_8u_C3R( src , srcStep , dst , dstStep , roiSize );
    __timer(&t1, &t2, "conversion"); // end clock

    // print dst array
    for(int c = 0; c < nChannels; c++){
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                printf("%d, ", dst[c*3*3 + i*3 + j]);
            }
        }
        printf("\n");
    }

    return 0; 
}
