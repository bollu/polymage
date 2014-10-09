#include <stdio.h>
#include <iostream>

#include <ippcore.h>
#include <ippvm.h>
#include <ipps.h>
#include <ippi.h>
#include <ippcv.h>

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "timing.h"

#include "vec_intrinsic.h"

#ifndef NRUNS
#define NRUNS 5
#endif

#ifndef SIZE
#define SIZE 512
#endif

using namespace cv;

void refHaar(int cols, int rows, float *input, float *&refDst_x_lo, float *&refDst_x_hi)
{
    for(int i = 0; i < rows; i++){
        #pragma ivdep
        for(int j = 0; j < cols/2; j++){
            refDst_x_lo[i*cols/2 + j] = (input[i*cols + 2*j  ] + input[i*cols + 2*j+1]) / 2;
            refDst_x_hi[i*cols/2 + j] = (input[i*cols + 2*j+1] - input[i*cols + 2*j  ]) / 2;
        }
    }
}

void refHaar_intr(int cols, int rows, float *input, float *&refDst_x_lo, float *&refDst_x_hi)
{
    /* load into vector reg and shuffle */

    const int imm_even = _MM_SHUFFLE(2,0,2,0);
    const int imm_odd = _MM_SHUFFLE(3,1,3,1);

    float half = 0.5f;
    VECTOR v_half;
    VECSPLAT(v_half, &half)
    
    int ubh = ((int)(cols/2)/VLEN)*(VLEN);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < ubh; j += VLEN){
            VECTOR first, second;
            VECTOR even, odd;
            VECTOR lo, hi;

            float * even_points = &input[i*cols + 2*j];
            float * odd_points  = &input[i*cols + 2*j + VLEN];

            first = VECLOADU(even_points)
            second = VECLOADU(odd_points)

            even = VECSHUFFLE(first, second, imm_even)
            odd = VECSHUFFLE(first, second, imm_odd)

            lo = VECADD(even, odd)
            hi = VECSUB(odd, even)

            lo = VECMUL(lo, v_half)
            hi = VECMUL(hi, v_half)

            VECSTOREU(lo, &refDst_x_lo[i*cols/2 + j])
            VECSTOREU(hi, &refDst_x_hi[i*cols/2 + j])
        }

        // epilogue
        for(int j = ubh; j < cols/2; j++){
            refDst_x_lo[i*cols/2 + j] = (input[i*cols + 2*j  ] + input[i*cols + 2*j+1]) / 2;
            refDst_x_hi[i*cols/2 + j] = (input[i*cols + 2*j+1] - input[i*cols + 2*j  ]) / 2;
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
    Mat imgIn;
    Rect roi(2000, 2000, SIZE, SIZE);
    if(argc == 2){
        imgIn = imread(argv[1]);
    }
    else{
        printf("Usage: %s image\n", argv[0]);
        exit(1);
    }
    image_info(imgIn);

    // cut the ROI
    Mat imgRoi = imgIn(roi).clone();
    int rows = imgRoi.size().height;
    int cols = imgRoi.size().width;

    // to grayscale
    Mat imgGray;
    cvtColor(imgRoi, imgGray, CV_BGR2GRAY);

    // to float
    Mat img;
    imgGray.convertTo(img, CV_32F, 1.0f/255);
    image_info(img);

    // result
    Mat dstCV(rows, cols, CV_32F);
    image_info(dstCV);

    // ===========================================================================================================

    // USING IPP

    // Initialize
    Mat img_ipp = img.clone();

    Mat dst_ipp_x_lo(rows, cols/2, CV_32F);
    Mat dst_ipp_x_hi(rows, cols/2, CV_32F);

	IppStatus status = ippStsNoErr;

    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    ippsWTHaarFwd_32f((const Ipp32f *)&img_ipp.data[0], rows*cols,
                      (Ipp32f *)&dst_ipp_x_lo.data[0],
                      (Ipp32f *)&dst_ipp_x_hi.data[0]);
    __timer(&t1, &t2, "ipp_x"); // end clock
    }
    printf("\n");

    // ===========================================================================================================

    // USING OPENCV
/*
    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    __timer(&t1, &t2, "ocv_x"); // end clock
    }
    printf("\n");
*/

    // ===========================================================================================================

    // REFERENCE

    float *refImg = (float *)malloc(sizeof(float) * rows * cols);
    float *refDst_x_lo = (float *)malloc(sizeof(float) * rows * cols/2);
    float *refDst_x_hi = (float *)malloc(sizeof(float) * rows * cols/2);

    // Initialize
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            refImg[i*cols + j] = img.at<float>(i, j);
        }
    }

    // Run
    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    refHaar(cols, rows, refImg, refDst_x_lo, refDst_x_hi);
    //refHaar_intr(cols, rows, refImg, refDst_x_lo, refDst_x_hi);
    __timer(&t1, &t2, "ref_x"); // end clock
    }
    printf("\n");

    // Collect
    Mat dstShow_ref(rows, cols, CV_32FC1);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols/2; j++){
            dstShow_ref.at<float>(i, j)        = refDst_x_lo[i*cols/2 + j];
            dstShow_ref.at<float>(i, j+cols/2) = refDst_x_hi[i*cols/2 + j];
        }
    }

    // ===========================================================================================================

    // DISPLAY
    #ifdef SHOW

    Mat dstShow_ipp(rows, cols, CV_32FC1);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols/2; j++)
            dstShow_ipp.at<float>(i, j) = dst_ipp_x_lo.at<float>(i, j);
        for(int j = 0; j < cols/2; j++)
            dstShow_ipp.at<float>(i, j+cols/2) = dst_ipp_x_hi.at<float>(i, j);
    }
    namedWindow( "IPP Result", WINDOW_NORMAL );

    //namedWindow( "CV Result", WINDOW_NORMAL );

    namedWindow( "Ref Result", WINDOW_NORMAL );

    for(;;) {
        int c;
        c = waitKey(10);
        if( (char)c == 27 )
        { break; }
        imshow( "IPP Result", dstShow_ipp );
        //imshow( "CV Result", imgShow_CV );
        imshow( "Ref Result", dstShow_ref );
    }
    #endif

    return 0; 
}
