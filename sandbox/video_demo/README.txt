Compiling the C++ file into a shared library
============================================

To compile all demo apps, run

$ make

To play the video demo, run the python script as

$ python video_demo.py path/to/video/file

The implementations work for a generic resolution.

For a sample video, try: https://peach.blender.org/download/ or
http://www.divx.com/en/devices/profiles/video


Options (key strokes)

h - to switch to harris corner detection
u - to switch to unsharp mask
b - to switch to bilateral grid
l - to switch to local laplacian filters

space - to toggle between OpenCV and PolyMage mode
        (applies only to 'harris')

n - to toggle between naive and optimized PolyMage code
