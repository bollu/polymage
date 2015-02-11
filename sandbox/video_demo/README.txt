Compiling the cpp file as a shared library
==========================================
icpc -xhost -openmp -fPIC -shared -o <file>.so <file>.cpp

icpc -xhost -openmp -fPIC -shared -o harris.so harris_polymage.cpp
icpc -xhost -openmp -fPIC -shared -o bilateral.so bilateral_polymage.cpp
icpc -xhost -openmp -fPIC -shared -o unsharp.so unsharp_polymage.cpp

Run the python script as

$ python video_demo.py path/to/video/file

The implementations are generated for video streams of resolution
1920x1080. Please use video files of the same resolution.

For a sample video try https://peach.blender.org/download/ or
http://www.divx.com/en/devices/profiles/video

Options in the demo
h - to switch to harris mode
b - to switch to bilateral mode
u - to switch to unsharp mode

space - to toggle OpenCV/PolyMage mode
        (only 'harris' for now)
