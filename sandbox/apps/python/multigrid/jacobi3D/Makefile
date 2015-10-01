CXX=icpc
CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared -DTIME
#CXX=g++
#CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize -fPIC -shared

all: new

new:
	python3.4 main.py 'new' 'W' 10

existing:
	python3.4 main.py 'existing' 'W' 10

ready:
	python3.4 main.py 'ready' 'W' 10

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *cycle.cpp
