#include "common.h"

using namespace cv;

void init(int argc, char *argv[], int &nruns, Mat &img1, Mat &img2, int &Rows, int &Cols)
{
    // Insufficient args
    if(argc < 3){
        std::cout <<"Usage:" <<std::endl
             <<argv[0] <<" image1 image2 [[#runs] [[height] [width]]]" <<std::endl;
        exit(1);
    }
    // > a.out image1 image2 [[#runs] [[height] [width]]]
    else{
        img1 = imread(argv[1]);
        img2 = imread(argv[2]);

        if(img1.rows != img2.rows || img1.cols != img2.cols){
            std::cout <<"Image Dimensions do not match!" <<std::endl;
            exit(1);
        }

        // optional arguments :
        // default
        nruns = 0;
        Rows = img1.rows;
        Cols = img1.cols;
        initOptionals(argc, argv, 3, nruns, Rows, Cols);
    }

    // Restrict the dimensions to the next lower power of 2
    Rows = (Rows & (Rows-1)) ? pow(2, ceil(log(Rows)/log(2))-1) : Rows;
    Cols = (Cols & (Cols-1)) ? pow(2, ceil(log(Cols)/log(2))-1) : Cols;
}

int main(int argc, char** argv)
{
    std::cout.precision(6);

    Mat img1, img2, res, img1f, img2f;
    int Rows = 0, Cols = 0, nruns;
    // Number of pyramid levels
    int L = 4;

    init(argc, argv, nruns, img1, img2, Rows, Cols);

    // Convert to floating point
    img1.convertTo(img1f, CV_32F, 1.0/255);
    img2.convertTo(img2f, CV_32F, 1.0/255);

    // Extract the region of interest
    Rect roi1(img1.cols/2-Cols/2, 0, Cols, Rows);
    Rect roi2(img2.cols/2-Cols/2, 0, Cols, Rows);
    Mat img1_region = img1f(roi1).clone();
    Mat img2_region = img2f(roi2).clone();

    image_info(img1_region);
    image_info(img2_region);

    // Mask to specify where to blend the images
    Mat mask(Rows, Cols, CV_32FC1);
    for (int i = 0; i < Rows; i++)
        for (int j = 0; j < Cols; j++)
            mask.at<float>(i, j) =  (j < (Cols/2)) ? 1.0 : 0.0;

    // Intermediate and final images
    Mat Gauss[2][L], Lapl[2][L], upsampled, 
        resLapl[L], maskGauss[L];

    for(int i = 0; i < nruns; i++) {
        Gauss[0][0] = img1_region;
        Gauss[1][0] = img2_region;
        maskGauss[0] = mask;
 
        TIMER__
        // Construct the gaussian pyramids of the source
        // images and the mask
        for (int l = 1; l < L; l++) {
            pyrDown(Gauss[0][l-1], Gauss[0][l]);
            pyrDown(Gauss[1][l-1], Gauss[1][l]);
            pyrDown(maskGauss[l-1], maskGauss[l]);
        }

        // Construct the laplacian pyramid of the source
        // images 
        for (int l = 0; l < L-1; l++) {
            pyrUp(Gauss[0][l+1], upsampled);
            Lapl[0][l] = Gauss[0][l] - upsampled;
            pyrUp(Gauss[1][l+1], upsampled);
            Lapl[1][l] = Gauss[1][l] - upsampled;
        }

        Lapl[0][L-1] = Gauss[0][L-1];
        Lapl[1][L-1] = Gauss[1][L-1];

        // Blend the laplacian pyramids using the mask 
        for (int l = 0; l < L; l++) {
            resLapl[l].create(Lapl[0][l].rows, Lapl[0][l].cols, Lapl[0][l].type());
            int rows = resLapl[l].rows;
            int cols = resLapl[l].cols;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    for (int c = 0; c < 3; c++) {
                        resLapl[l].at<Vec3f>(i, j)[c] = 
                               Lapl[0][l].at<Vec3f>(i, j)[c] * maskGauss[l].at<float>(i, j)
                             + Lapl[1][l].at<Vec3f>(i, j)[c] * (1 - maskGauss[l].at<float>(i, j));
                    }
                }
            }
        }

        // Collapse the pyramid
        for(int l = L-2; l >=0; l--) {
            pyrUp(resLapl[l+1], upsampled);
            resLapl[l] = resLapl[l] + upsampled;
        }
        __TIMER("OpenCV")
    }

    #ifdef SHOW
    // Create window handles
    const char* window_image1 = "Image 1";
    const char* window_image2 = "Image 2";
    const char* window_opencv_result = "OpenCV Result";
    // Create windows
    namedWindow( window_image1, WINDOW_NORMAL );
    namedWindow( window_image2, WINDOW_NORMAL );
    namedWindow( window_opencv_result, WINDOW_NORMAL );
    // Display
    for(;;) {
        int c;
        c = waitKey(10);
        if( (char)c == 27 )
        { break; }
        imshow( window_image1, img1_region);
        imshow( window_image2, img2_region);
        imshow( window_opencv_result, resLapl[0]);
    }
    #endif

    return 0;
}
