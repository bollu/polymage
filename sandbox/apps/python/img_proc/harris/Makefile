#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared -DTIME
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize -fPIC -shared

all: new

new:
	python3.4 main.py --mode='new' --img='../../../../../images/venice_wikimedia.jpg' --graph-gen --runs=1 -t --display

tune:
	python3.4 main.py --mode='tune' --img='../../../../../images/venice_wikimedia.jpg' --runs=1 -t --display

existing:
	python3.4 main.py --mode='existing' --img='../../../../../images/venice_wikimedia.jpg' --runs=1 -t --display

ready:
	python3.4 main.py --mode='ready' --img='../../../../../images/venice_wikimedia.jpg' --runs=1 -t --display

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *.cpp
