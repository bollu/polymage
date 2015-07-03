CXX=icpc
CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared -DTIME
#CXX=g++
#CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize -fPIC -shared

all: new

new:
	python main.py 'new' 10

existing:
	python main.py 'existing' 10

ready:
	python main.py 'ready' 10

clean:
	rm -f *.pyc *.so *.dot *.png

cleaner:
	rm -f *.pyc *.so *.dot *.png vcycle.cpp
