#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace cv;

#ifdef TIME
    double t_start, t_end;
    #define TIMER__ t_start = getTickCount();
    #define __TIMER(section) t_end = getTickCount(); \
                    std::cout <<section <<" time = " \
                    <<1000 * ((t_end - t_start)/getTickFrequency()) <<" ms" \
                    <<std::endl;
#else
    #define TIMER__
    #define __TIMER(section)
#endif

inline void image_info(Mat &img)
{
    std::cout << "Channels: " << img.channels() << std::endl;
    std::cout << "Depth: "  << img.type() << std::endl;     
    std::cout << "Size: "  << img.size() << std::endl;     
}

inline void setmin(int &a, const int b)
{
    a = (a > b) ? a : b;
}

inline void setmax(int &a, const int b)
{
    a = (a > b) ? b : a;
}

void clamp(int &a, const int b, const int c)
{
    a = (a < b) ? b : (a > c) ? c : a;
}

void initOptionals(const int argc, char *argv[], const int offset, int &nruns, int &Rows, int &Cols)
{
    int tempRows = Rows, tempCols = Cols;

    if(argc > offset)
        nruns = atoi(argv[offset]);
    if(argc > offset+1){
        // height of roi image
        tempRows = atoi(argv[offset+1]);
        if(argc > offset+2)
            // width of roi image
            tempCols = atoi(argv[offset+2]);
        else
            tempCols = tempRows;
    }
    // Clamp the parameters
    setmin(nruns, 1);
    clamp(tempRows, 1, Rows);
    clamp(tempCols, 1, Cols);

    Rows = tempRows;
    Cols = tempCols;
}
