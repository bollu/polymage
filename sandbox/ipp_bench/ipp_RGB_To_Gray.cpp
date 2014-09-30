#include <stdio.h>

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

	Ipp8u src[12*3] = {
		255, 0, 0, 255, 0, 0, 255, 0, 0, 255, 0, 0, 
		0, 255, 0, 0, 255, 0, 0, 255, 0, 0, 255, 0,
		0, 0, 255, 0, 0, 255, 0, 0, 255, 0, 0, 255};

	Ipp8u dst[4*3];
	IppiSize srcRoi = { 4, 3 };

    timer__(&t1); // start clock

    // conversion from RGB to Gray
	ippiRGBToGray_8u_C3C1R ( src, 12, dst, 4, srcRoi );

    __timer(&t1, &t2, "conversion"); // end clock

    // print dst array
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 4; j++){
            printf("%d, ", dst[i*4 + j]);
        }
        printf("\n");
    }

    return 0; 
}
