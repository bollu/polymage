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

void refHist(int cols, int rows, unsigned char *refImg, int *& refDst)
{
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            refDst[refImg[i*cols + j]]++;
        }
    }
}

extern void refHarris(int cols, int rows, unsigned char *input, unsigned char *& harris);

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

    // to CV_8U
    Mat img8u;
    imgRoi.convertTo(img8u, CV_8U, 1);

    // to grayscale
    Mat img;
    cvtColor(img8u, img, CV_BGR2GRAY);

    // Histogram Properties
    const int nBins = 256;
    int nLevels[] = { nBins };


    // ===========================================================================================================
    // USING IPP

    // Initialize

    Ipp32f lowerLevel[] = {0};
    Ipp32f upperLevel[] = {256};
    Ipp32f pLevels[nBins], *ppLevels[1];

    int sizeHistObj, sizeBuffer;
    IppiHistogramSpec* pHistObj;
    Ipp8u* pBuffer;
    Ipp32u pHistVec[nBins-1];

    IppiSize ippRoi = {rows, cols};

    ippiHistogramGetBufferSize(ipp8u, ippRoi, nLevels, 1/*nChannels*/, 1/*Uniform Bins*/, &sizeHistObj, &sizeBuffer);
    pHistObj = (IppiHistogramSpec*)ippsMalloc_8u( sizeHistObj );
    pBuffer = (Ipp8u*)ippsMalloc_8u( sizeBuffer );

    ippiHistogramUniformInit( ipp8u, lowerLevel, upperLevel, nLevels, 1/*nChannels*/, pHistObj );

    ppLevels[0] = pLevels;
    IppStatus sts = ippiHistogramGetLevels( pHistObj, ppLevels );

    // Calculate

    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    sts = ippiHistogram_8u_C1R( (const Ipp8u *)&img.data[0], cols, ippRoi, pHistVec, pHistObj, pBuffer );
    __timer(&t1, &t2, "ipp"); // end clock
    }
    printf("\n");

    // Free

    ippsFree( pHistObj );
    ippsFree( pBuffer );


    // ===========================================================================================================
    // USING OPENCV

    // Initialize
    Mat dstCV;
    float range[2] = {0, nBins};
    const float *histRange = { range };

    // Compute
    for(int runs = 0; runs < NRUNS; runs++){
    timer__(&t1); // start clock
    calcHist(&img, 1/*nImages*/, 0/*channels*/, Mat()/*Mask*/, dstCV, 1, &nBins, &histRange, true/*uniform*/, false/*accumulate*/);
    __timer(&t1, &t2, "opencv"); // end clock
    }
    printf("\n");


    // ===========================================================================================================
    // REFERENCE

    // Initialize
    unsigned char *refImg = (unsigned char *)malloc(sizeof(unsigned char) * (rows) * (cols));
    int *refDst = (int *)malloc(sizeof(int) * nBins);

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            refImg[i*cols + j] = img.at<unsigned char>(i, j);
        }
    }

    for(int runs = 0; runs < NRUNS; runs++){
    for(int b = 0; b < nBins; b++)
        refDst[b] = 0;

    // Compute
    timer__(&t1); // start clock
    refHist(cols, rows, refImg, refDst);
    __timer(&t1, &t2, "ref"); // end clock

    }
    printf("\n");

    // ===========================================================================================================
    // DISPLAY

    #ifdef SHOW

    // Draw the histograms

    int hist_w = 512; int hist_h = 400;
    int bin_w = cvRound( (double) hist_w/nBins );

    // ipp Mat
    Mat dstIpp(1, nBins, CV_32F);
    for(int b = 0; b < nBins; b++)
        dstIpp.at<float>(0, b) = pHistVec[b];

    // ref Mat
    Mat dstRef(1, nBins, CV_32F);
    for(int b = 0; b < nBins; b++)
        dstRef.at<float>(0, b) = refDst[b];

    /*
    int count = 0;
    for(int b = 0; b < nBins; b++){
        printf("%d, %f, %d\n", dstIpp.at<int>(0, b), dstCV.at<float>(0, b), dstRef.at<int>(0, b));
        count += dstIpp.at<int>(0, b);
    }
    printf("\ntotal: %d\n", count);
    */

    // Create the display Mat
    Mat histShowIPP( hist_h, hist_w, CV_8UC3, Scalar( 0, 0, 0 ) );
    Mat histShowCV ( hist_h, hist_w, CV_8UC3, Scalar( 0, 0, 0 ) );
    Mat histShowRef( hist_h, hist_w, CV_8UC3, Scalar( 0, 0, 0 ) );

    // Normalize the result to [ 0, histImage.rows ]
    normalize(dstIpp, dstIpp, 0, histShowIPP.rows, NORM_MINMAX, -1, Mat() );
    normalize(dstCV , dstCV , 0, histShowCV.rows , NORM_MINMAX, -1, Mat() );
    normalize(dstRef, dstRef, 0, histShowRef.rows, NORM_MINMAX, -1, Mat() );

    // Draw IPP
    for( int i = 1; i < nBins; i++ )
    {
        line( histShowIPP, Point( bin_w*(i-1), hist_h - cvRound(dstIpp.at<float>(i-1)) ) ,
                           Point( bin_w*(i)  , hist_h - cvRound(dstIpp.at<float>(i)) ),
                           Scalar( 0, 0, 255 ), 2, 8, 0  );
    }

    // Draw CV
    for( int i = 1; i < nBins; i++ )
    {
        line( histShowCV , Point( bin_w*(i-1), hist_h - cvRound(dstCV.at<float>(i-1)) ) ,
                           Point( bin_w*(i)  , hist_h - cvRound(dstCV.at<float>(i)) ),
                           Scalar( 0, 0, 255 ), 2, 8, 0  );
    }

    // Draw REF
    for( int i = 1; i < nBins; i++ )
    {
        line( histShowRef, Point( bin_w*(i-1), hist_h - cvRound(dstRef.at<float>(i-1)) ) ,
                           Point( bin_w*(i)  , hist_h - cvRound(dstRef.at<float>(i)) ),
                           Scalar( 0, 0, 255 ), 2, 8, 0  );
    }


    namedWindow( "IPP Result", WINDOW_NORMAL );
    namedWindow( "CV Result", WINDOW_NORMAL );
    namedWindow( "Ref Result", WINDOW_NORMAL );

    for(;;) {
        int c;
        c = waitKey(10);
        if( (char)c == 27 )
        { break; }
        imshow( "IPP Result", histShowIPP );
        imshow( "CV Result", histShowCV );
        imshow( "Ref Result", histShowRef );
    }
    #endif

    return 0; 
}
