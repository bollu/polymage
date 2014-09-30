(*) set ipp's library to LD_LIBRARY_PATH
(*) make note of deprecated routines as and when encountered in our target
(*) template compile string :
    icpc -O3 -xhost -openmp
    -I [path-to-intel-studio]/intel/studio/ipp/include]
    -lippcore -lippvm -lipps -lippi -lippcc
