#include <stdio.h>
#include <iostream>

#include <ippcore.h>
#include <ippvm.h>
#include <ipps.h>
#include <ippi.h>
#include <ippcc.h>

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "timing.h"

#ifndef NRUNS
#define NRUNS 5
#endif

#ifndef SIZE
#define SIZE 512
#endif

using namespace cv;

void refGray(int cols, int rows, const float *img, float *& res)
{
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < rows; i++) {
        #pragma simd
        for(int j = 0; j < cols; j++) {
            res[i*cols + j] = (img[i*cols*3 + j*3 + 0] * 0.114f +
                               img[i*cols*3 + j*3 + 1] * 0.587f +
                               img[i*cols*3 + j*3 + 2] * 0.299f);
            
            /*res[i*cols + j] = img[0*rows*cols + i*cols + j] * 0.114f +
                              img[1*rows*cols + i*cols + j] * 0.587f +
                              img[2*rows*cols + i*cols + j] * 0.299f;*/
        }
    }
}

void image_info(Mat &img)
{
    std::cout << "Channels: " << img.channels() << std::endl;
    std::cout << "Depth: "  << img.type() << std::endl;
    std::cout << "Size: "  << img.size() << std::endl;
}

int main(int argc, char *argv[])
{
    cvUseOptimized(1);
    timeval t1, t2;

    // FIXME (just for attention - set path and roi properly)
    Mat imgf;
    Rect roi(0, 0, SIZE, SIZE);
    if(argc == 2){
        imgf = imread(argv[1]);
    }
    else{
        printf("Usage: %s image\n", argv[0]);
        exit(1);
    }

    // to float
    Mat img1;
    imgf.convertTo(img1, CV_32F, 1.0f/255);
    Mat img = img1(roi).clone();
    Mat dst(img.size().height, img.size().width, CV_32F);

    // size
    int rows = img.size().height;
    int cols = img.size().width;

    image_info(img);
    image_info(dst);

	IppiSize imgRoi = { dst.size().width, dst.size().height };

    // ===========================================================================================================

    // USING IPP

    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
	ippiRGBToGray_32f_C3C1R ( (const Ipp32f *)&img.data[0], 3*img.size().width, (Ipp32f *)&dst.data[0], dst.size().width, imgRoi );
    __timer(&t1, &t2, "ipp"); // end clock
    }
    printf("\n");

    // ===========================================================================================================

    // USING OPENCV

    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    cv::cvtColor(img, dst, CV_BGR2GRAY);
    __timer(&t1, &t2, "ocv"); // end clock
    }
    printf("\n");

    // ===========================================================================================================

    // REFERENCE

    // Initialize
    float *refImg = (float *)malloc(sizeof(float) * 3 * rows * cols);
    float *refDst = (float *)malloc(sizeof(float) * rows * cols);

    for(int c = 0; c < 3; c++)
        //#pragma omp parallel for
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                refImg[i*cols*3 + j*3 + c] = img.at<Vec3f>(i, j)[c];
                //refImg[c*cols*rows + i*cols + j] = img.at<Vec3f>(i, j)[c];

    // Run
    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    refGray(cols, rows, refImg, refDst);
    __timer(&t1, &t2, "ref"); // end clock
    }
    printf("\n");

    // Collect
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            dst.at<float>(i, j) = refDst[i*cols + j];

    // ===========================================================================================================

    // DISPLAY

    #ifdef SHOW
    namedWindow( "Image", WINDOW_NORMAL );
    namedWindow( "Result", WINDOW_NORMAL );
    for(;;) {
        int c;
        c = waitKey(10);
        if( (char)c == 27 )
        { break; }
        imshow( "Image", img );
        imshow( "Result", dst );
    }
    #endif

    return 0; 
}
