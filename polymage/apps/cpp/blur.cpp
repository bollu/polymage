#include "common.h"

using namespace cv;

void init(int argc, char *argv[], int &nruns, Mat &img, int &Rows, int &Cols)
{
    // Insufficient args
    if(argc < 2){
        std::cout <<"Usage:" <<std::endl
             <<argv[0] <<" image [[#runs] [[height] [width]]]" <<std::endl;
        exit(1);
    }
    // > a.out image [[#runs] [[height] [width]]]
    else{
        img = imread(argv[1]);

        // optional arguments :
        // default
        nruns = 0;
        Rows = img.rows;
        Cols = img.cols;
        initOptionals(argc, argv, 2, nruns, Rows, Cols);
    }
}

int main(int argc, char** argv)
{    
    std::cout.precision(6);

    // Input and intermediate matrices
    Mat img, imgf;
    Mat blurx, blury;
    int Rows = 0, Cols = 0, nruns;

    init(argc, argv, nruns, img, Rows, Cols);

    // Convert to floating point
    img.convertTo(imgf, CV_32F, 1.0/255);

    // Extract region of interest
    Rect roi(0, 0, Cols, Rows);
    Mat img_region = imgf(roi).clone();

    image_info(img_region);

    // Blur kernels
    Mat kernx = (Mat_<float>(1,3) <<  1.0/3, 1.0/3, 1.0/3);
    Mat kerny = (Mat_<float>(3,1) <<  1.0/3, 1.0/3, 1.0/3);

    // Compute
    for (int i = 0; i < nruns; i++) {
        TIMER__
        filter2D(img_region, blurx, img_region.depth(), kernx, Point(-1,-1), BORDER_CONSTANT); 
        filter2D(blurx, blury, blurx.depth(), kerny, Point(-1,-1), BORDER_CONSTANT);  
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
        imshow( window_opencv_result, blury);
    }
    #endif

    return 0;
}
