CXX=icpc
CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared -DTIME
#CXX=g++
#CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize -fPIC -shared

all: new

new:
	python3.4 main.py 'new' 'V' 10

tune:
	python3.4 main.py 'tune' 'V' 10

existing:
	python3.4 main.py 'existing' 'V' 10

ready:
	python3.4 main.py 'ready' 'V' 10

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *cycle.cpp
