CXX=icpc
CXX_FLAGS=-openmp -O3 -xhost -fPIC -shared
#CXX=g++
#CXX_FLAGS=-fopenmp -O3 -march=native -fPIC -shared

all: norm.so new

new:
	python main.py 'D' 'new'

existing:
	python main.py 'D' 'existing'

ready:
	python main.py 'D' 'ready'

norm.so: norm.cpp
	${CXX} ${CXX_FLAGS} norm.cpp -o norm.so

clean:
	rm -f *.pyc *.so *.out *.png

cleaner:
	rm -f *.pyc *.so *.out *.png *.dot resid.cpp mgU_mgR_.cpp
