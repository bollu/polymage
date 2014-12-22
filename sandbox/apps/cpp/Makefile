CXX=icpc
CXX_FLAGS=-ansi-alias -ipo -openmp -O3 -xhost
LIBS=-L /usr/local/lib/ 
LD_FLAGS=-lopencv_imgproc -lopencv_core -lopencv_highgui
APP_FLAGS=-DTIME -DSHOW

all: blur harris pyramid_blend unsharp_mask

blur: blur.o
	./blur ~/images/room1.png 5 1920 1080

harris: harris.o
	./harris ~/images/port.jpg 75 5 6400 6400

pyramid_blend: pyramid_blend.o
	./pyramid_blend ~/images/room1.png ~/images/room2.png 5 1024 512

unsharp_mask: unsharp_mask.o
	./unsharp_mask ~/images/chip.jpg 0.001 3 5 5000 4000

blur.o:
	$(CXX) $(CXX_FLAGS) $(LIBS) $(LD_FLAGS) blur.cpp -o blur $(APP_FLAGS)

harris.o:
	$(CXX) $(CXX_FLAGS) $(LIBS) $(LD_FLAGS) harris.cpp -o harris $(APP_FLAGS)

pyramid_blend.o:
	$(CXX) $(CXX_FLAGS) $(LIBS) $(LD_FLAGS) pyramid_blending.cpp -o pyramid_blend $(APP_FLAGS)

unsharp_mask.o:
	$(CXX) $(CXX_FLAGS) $(LIBS) $(LD_FLAGS) unsharp_mask.cpp -o unsharp_mask $(APP_FLAGS)

clean:
	rm -f blur blur.o harris harris.o pyramid_blend pyramid_blend.o unsharp_mask unsharp_mask.o
