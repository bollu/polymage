#include "common.h"

using namespace cv;

void init(int argc, char *argv[], int &nruns, Mat &img, float &threshold, int &Rows, int &Cols)
{
    // Insufficient args
    if(argc < 3){
        std::cout <<"Usage:" <<std::endl
             <<argv[0] <<" image threshold [[#runs] [[height] [width]]]" <<std::endl;
        exit(1);
    }
    // > a.out image threshold [[#runs] [[height] [width]]]
    else{
        img = imread(argv[1]);
        threshold = atof(argv[2]);

        // optional arguments :
        // default
        nruns = 0;
        Rows = img.rows;
        Cols = img.cols;
        initOptionals(argc, argv, 3, nruns, Rows, Cols);
    }
}

int main(int argc, char** argv)
{    
    std::cout.precision(6);

    // Input, Output and Intermediate matrices
    Mat img, imgGray, imgGrayf;
    Mat opencv_harris;
    int Rows, Cols, nruns;

    // corner threshold
    float threshold = 0.0f;

    init(argc, argv, nruns, img, threshold, Rows, Cols);

    // Convert to grayscale
    cv::cvtColor(img, imgGray, CV_BGR2GRAY);

    // Convert to floating point
    imgGray.convertTo(imgGrayf, CV_32F, 1.0/255);

    // Extract region of interest
    Rect roi(0, 0, Cols, Rows);
    Mat img_region = imgGrayf(roi).clone();

    image_info(img_region);

    // Image base to display
    Mat img_show = img(roi).clone();
    Mat img_show_harris = img(roi).clone();

    // Compute corner map
    for (int i = 0; i < nruns; i++) {
    	TIMER__	
        cornerHarris(img_region, opencv_harris, 3, 3, 0.04);
	    __TIMER("OpenCV")
    }
    normalize(opencv_harris, opencv_harris, 0, 255, NORM_MINMAX, CV_32FC1, Mat() );
    convertScaleAbs(opencv_harris, opencv_harris);

	#ifdef SHOW
    //Create window handles
    const char* window_image = "Image";
    const char* window_opencv_result = "OpenCV Result";
    // Encircle the pixels at which intensity exceeds corner threshold
    for( int i = 0; i < Rows; i++ ){ 
        for( int j = 0; j < Cols; j++ ){
            if( opencv_harris.at<unsigned char>(i, j) > threshold)
                circle( img_show_harris, Point( j, i ), 5,  Scalar(0, 0, 255), 2, 8, 0 );
        }
    }
    // Create windows
    namedWindow( window_image, WINDOW_NORMAL );
    namedWindow( window_opencv_result, WINDOW_NORMAL );
    // Display
    for(;;) {
        int c;
        c = waitKey(10);
        if( (char)c == 27 )
        { break; }
        imshow( window_image, img_show);
        imshow( window_opencv_result, img_show_harris);
    }
	#endif

    return 0;
}
