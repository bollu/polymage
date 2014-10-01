(*) set ipp's library to LD_LIBRARY_PATH
(*) make note of deprecated routines as and when encountered in our target
(*) template compile string :
    icpc -O3 -xhost -openmp
    -I [path-to-intel-studio]/intel/studio/ipp/include]
    -L /usr/local/lib/ 
    -lippcore -lippvm -lipps -lippi -lippcc
    -lopencv_imgproc -lopencv_core -lopencv_highgui

