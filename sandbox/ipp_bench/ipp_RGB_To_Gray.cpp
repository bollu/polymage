// compiler says "ippiRGBToGray_8u_C3C1R" undefined keyword

#include <stdio.h>

#include <ippcore.h>
#include <ippvm.h>
#include <ipps.h>
#include <ippi.h>

# define nChannels 3

int main ()
{
	Ipp8u src[12*3] = {
		255, 0, 0, 255, 0, 0, 255, 0, 0, 255, 0, 0, 
		0, 255, 0, 0, 255, 0, 0, 255, 0, 0, 255, 0,
		0, 0, 255, 0, 0, 255, 0, 0, 255, 0, 0, 255};

	Ipp8u dst[4*3];
	IppiSize srcRoi = { 4, 3 };

	ippiRGBToGray_8u_C3C1R ( src, 12, dst, 4, srcRoi );

    return 0; 
}
