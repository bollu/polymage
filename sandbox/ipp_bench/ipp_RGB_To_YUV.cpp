// compiler says "ippiRGBToYUV_8u_C3R" undefined keyword

#include <stdio.h>

#include <ippcore.h>
#include <ippvm.h>
#include <ipps.h>
#include <ippi.h>

# define nChannels 3

int main ()
{
    Ipp8u src [3*3*nChannels] = {
        255, 0, 0, 255, 0, 0, 255, 0, 0,
        0, 255, 0, 0, 255, 0, 0, 255, 0,
        0, 0, 255, 0, 0, 255, 0, 0, 255};

    Ipp8u dst [3*3* nChannels ];

    IppiSize roiSize = { 3, 3 };
    IppStatus st = ippStsNoErr ;

    int srcStep = 3 * nChannels ;
    int dstStep = 3 * nChannels ; 

    st = ippiRGBToYUV_8u_C3R( src , srcStep , dst , dstStep , roiSize );

    if ( st == ippStsNoErr)
        printf("\n ************* passed ****************\n");
    else
        printf("\n ************* failed ****************\t");

    return 0; 
}
