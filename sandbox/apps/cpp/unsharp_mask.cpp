#include "common.h"

using namespace cv;

void init(const int argc, char *argv[], int &nruns, Mat &img, double &threshold, double &weight, int &Rows, int &Cols)
{
    // Insufficient args
    if(argc < 4){
        std::cout <<"Usage:" <<std::endl
             <<argv[0] <<" image threshold weight [[#runs] [[height] [width]]]" <<std::endl;
        exit(1);
    }
    // > a.out image [[#runs] [[height] [width]]]
    else{
        img = imread(argv[1]);
        threshold = (double)atof(argv[2]);
        weight = (double)atof(argv[3]);

        // optional arguments :
        // default
        nruns = 0;
        Rows = img.rows;
        Cols = img.cols;
        initOptionals(argc, argv, 4, nruns, Rows, Cols);
    }
}

int main(int argc, char** argv)
{    
    std::cout.precision(6);

    // Input and intermediate matrices 
    Mat img, imgf;
    Mat blury, mask, ocv_sharpened;
    int Rows, Cols, nruns;
    double threshold, weight, sigma = -0.75;
    // threshold = 0.001, weight = 3;

    init(argc, argv, nruns, img, threshold, weight, Rows, Cols);

    // Convert to floating point
    img.convertTo(imgf, CV_32F, 1.0/255);

    // Extract region of interest
    Rect roi(0, 0, Cols, Rows);
    Mat img_region = imgf(roi).clone();

    image_info(img_region);

    // Compute the sharpened image
    for (int i = 0; i < nruns; i++) {
        TIMER__
        GaussianBlur(img_region, blury, Size(5, 5), sigma, sigma);
        mask = (abs(img_region - blury) < threshold);
        ocv_sharpened = img_region*(1+weight) + blury*(-weight);
        img_region.copyTo(ocv_sharpened, mask);
        __TIMER("OpenCV")
    }

    #ifdef SHOW
    // Create window handles
    const char* window_image  = "Image";
    const char* window_opencv_result  = "OpenCV Result";
    // Create windows
    namedWindow( window_image, WINDOW_NORMAL );
    namedWindow( window_opencv_result, WINDOW_NORMAL );
    // Display
    for(;;) {
        int c;
        c = waitKey(10);
        if( (char)c == 27 )
        { break; }
        imshow( window_image, img_region);
        imshow( window_opencv_result, ocv_sharpened);
    }
    #endif

    return 0;
}
