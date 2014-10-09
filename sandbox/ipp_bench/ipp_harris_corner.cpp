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

#ifndef NRUNS
#define NRUNS 5
#endif

#ifndef SIZE
#define SIZE 4096
#endif

using namespace cv;

extern void refHarris(int cols, int rows, float *input, float *& harris);

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
    Mat dst(rows, cols, CV_32F);
    image_info(dst);

    // ===========================================================================================================

    // USING IPP

    // Initialize
    Mat img_ipp(rows+2, cols+2, CV_32F, cvScalar(0.0f));
    img.copyTo(img_ipp(Rect(1, 1, cols, rows))); // source with ghost

    Mat dst_ipp = dst.clone();

	IppiSize imgRoi_ipp = { cols, rows };

    int bufSize = 0;
    Ipp8u* pBuffer = NULL;
	IppStatus status = ippStsNoErr;

    double scale = std::pow(12.0, -4.0);

	// Compute work buffer size
	status = ippiHarrisCornerGetBufferSize(imgRoi_ipp, ippMskSize3x3,          3, ipp32f,        1,  &bufSize);
	       //ippiHarrisCornerGetBufferSize(       ROI,   filter mask, avgWndSize,   type, channels, allocated);

	// Allocate Memory
	if (status == ippStsNoErr) pBuffer = ippsMalloc_8u(bufSize);

    // Compute
	if (pBuffer != NULL){
        for(int runs = 0; runs < NRUNS; runs++){
        timer__(&t1); // start clock
        status = ippiHarrisCorner_32f_C1R((const Ipp32f *)&(img_ipp.data[1*(cols+2)+1]), (cols+2)*sizeof(Ipp32f), (Ipp32f *)&dst_ipp.data[0], cols*sizeof(Ipp32f), imgRoi_ipp,
                                          ippFilterSobel, ippMskSize3x3, 3,
                                          0.04f, (Ipp32f)scale, ippBorderConst, 0.0f, pBuffer);
		//       ippiHarrisCorner_32f_C1R(pSrc, srcStep, pDst, dstStep, roiSize,
        //                                filterType, filterMask, avgWndSize,
        //                                k, scale, borderType, borderValue, pBuffer);
        __timer(&t1, &t2, "ipp"); // end clock
        std::cout <<"\tstatus: " <<ippGetStatusString(status) <<std::endl;
        }
        ippsFree(pBuffer);
    }
    printf("\n");

    // ===========================================================================================================

    // USING OPENCV

    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    cornerHarris(img, dst, 3, 3, 0.04);
    __timer(&t1, &t2, "opencv"); // end clock
    }
    printf("\n");

    // ===========================================================================================================

    // REFERENCE

    Mat dstRef = dst.clone();
    float *refImg = (float *)malloc(sizeof(float) * (rows+2) * (cols+2));
    float *refDst = (float *)malloc(sizeof(float) * (rows+2) * (cols+2));

    // Initialize
    //#pragma omp parallel for
    for(int i = 0; i < rows+2; i++){
        for(int j = 0; j < cols+2; j++){
            int ii = i - 1 > 0 ?  (i -1 < rows ? i-1:  rows-1): 0;
            int jj = j - 1 > 0 ?  (j -1 < cols ? j-1:  cols-1): 0;
            refImg[(i)*(cols+2) + (j)] = img.at<float>(ii, jj);
        }
    }

    // Run
    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    refHarris(cols, rows, refImg, refDst);
    __timer(&t1, &t2, "ref"); // end clock
    }
    printf("\n");

    // Collect
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            dstRef.at<float>(i, j) = refDst[(i+1)*(cols+2) + (j+1)];

    // ===========================================================================================================

    // DISPLAY
    #ifdef SHOW

    /*
    // IPP's
    Mat imgShow_IPP = imgRoi.clone();
    //normalize(dst_ipp, dst_ipp, 0, 255, NORM_MINMAX, CV_32FC1, Mat() );
    //convertScaleAbs(dst_ipp, dst_ipp);
    for( int i = 0; i < dst_ipp.rows ; i++ )
        for( int j = 0; j < dst_ipp.cols; j++ )
            if( dst_ipp.at<unsigned char>(i, j) > 0)
                circle( imgShow_IPP, Point( j, i ), 5,  Scalar(0, 0, 255), 2, 8, 0 );

    // OpenCV's
    Mat imgShow_CV = imgRoi.clone();
    //normalize(dst, dst, 0, 255, NORM_MINMAX, CV_32FC1, Mat() );
    //convertScaleAbs(dst, dst);
    for( int i = 0; i < dst.rows ; i++ )
        for( int j = 0; j < dst.cols; j++ )
            if( dst.at<unsigned char>(i, j) > 75)
                circle( imgShow_CV, Point( j, i ), 5,  Scalar(0, 0, 255), 2, 8, 0 );

    // Ref's
    Mat imgShow_Ref = imgRoi.clone();
    //normalize(dstRef, dstRef, 0, 255, NORM_MINMAX, CV_32FC1, Mat() );
    //convertScaleAbs(dstRef, dstRef);
    for( int i = 0; i < dstRef.rows; i++ )
        for( int j = 0; j < dstRef.cols; j++ )
            if( dstRef.at<unsigned char>(i, j) > 75)
                circle( imgShow_Ref , Point( j+1, i+1 ), 5,  Scalar(0, 0, 255), 2, 8, 0 );
    */

    namedWindow( "IPP Result", WINDOW_NORMAL );

    namedWindow( "CV Result", WINDOW_NORMAL );

    namedWindow( "Ref Result", WINDOW_NORMAL );

    for(;;) {
        int c;
        c = waitKey(10);
        if( (char)c == 27 )
        { break; }
        //imshow( "IPP Result", imgShow_IPP );
        imshow( "IPP Result", dst_ipp );
        //imshow( "CV Result", imgShow_CV );
        imshow( "CV Result", dst );
        //imshow( "Ref Result", imgShow_Ref );
        imshow( "Ref Result", dstRef );
    }
    #endif

    return 0; 
}
