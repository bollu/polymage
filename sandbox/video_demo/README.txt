Compiling the cpp file as a shared library
==========================================
icpc -xhost -openmp -fPIC -shared -o <file>.so <file>.cpp

icpc -xhost -openmp -fPIC -shared -o harris.so harris_polymage.cpp
icpc -xhost -openmp -fPIC -shared -o bilateral.so bilateral_polymage.cpp

Run the python script as

python video_benchmark.py path_to_video_file 

The bilateral grid  implementation is generated for a video stream of 
resolution 1920x1080; so please use a video file of that resolution. 

For a sample video try https://peach.blender.org/download/ or
http://www.divx.com/en/devices/profiles/video

Options in the demo
h - to switch to harris mode
space - to toggle OpenCV/PolyMage mode
b - to switch to bilateral mode
