Compiling the cpp file as a shared library
==========================================
To compile all the apps, run

$ make

To play the video demo, run the python script as

$ python video_demo.py path/to/video/file

The implementations are generated for generic resolution.

For a sample video try https://peach.blender.org/download/ or
http://www.divx.com/en/devices/profiles/video

Options in the demo
h - to switch to harris mode
b - to switch to bilateral mode
u - to switch to unsharp mode

space - to toggle OpenCV/PolyMage mode
        (only 'harris' for now)
n - to toggle Naive/Opt PolyMage code
