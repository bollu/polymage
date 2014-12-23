CXX=icpc
CXXFLAGS=-xhost -openmp -fPIC -shared

all: harris bilateral

harris:
	$(CXX) $(CXXFLAGS) -o harris.so harris_polymage.cpp

bilateral:
	$(CXX) $(CXXFLAGS) -o bilateral.so bilateral_polymage.cpp

clean:
	rm -f *.so
