#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared -DTIME
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize -fPIC -shared
PYTHON=python3

all: new

new:
	$(PYTHON) main.py --mode='new' --img='../../../../../images/venice_wikimedia.jpg' --graph-gen --runs=1 -t --display --rowdiff=2832 --coldiff=4256 --off_left=24 --total_pad=56

tune:
	$(PYTHON) main.py --mode='tune' --img='../../../../../images/venice_wikimedia.jpg' --runs=1 -t --display --rowdiff=2832 --coldiff=4256 --off_left=24 --total_pad=56

existing:
	$(PYTHON) main.py --mode='existing' --img='../../../../../images/venice_wikimedia.jpg' --runs=1 -t --display --rowdiff=2832 --coldiff=4256 --off_left=24 --total_pad=56

ready:
	$(PYTHON) main.py --mode='ready' --img='../../../../../images/venice_wikimedia.jpg' --runs=1 -t --display --rowdiff=2832 --coldiff=4256 --off_left=24 --total_pad=56

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *.cpp
