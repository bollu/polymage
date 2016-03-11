#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost -ansi-alias
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize
SHARED_LIB_FLAGS=-fPIC -shared

PYTHON=python3.4

CLASS='W'  # NAS-MG Problem Class
RUNS=1

PROB_ARGS=--class=${CLASS} --runs=${RUNS}

all: norm.so new

new:
	${PYTHON} main.py --mode='new' ${PROB_ARGS} --cxx="${CXX}" --cxx_flags="${CXX_FLAGS}" -t

tune:
	${PYTHON} main.py --mode='tune' ${PROB_ARGS} --cxx="${CXX}" --cxx_flags="${CXX_FLAGS}" -t

existing:
	${PYTHON} main.py --mode='existing' ${PROB_ARGS} --cxx="${CXX}" --cxx_flags="${CXX_FLAGS}" -t

ready:
	${PYTHON} main.py --mode='ready' ${PROB_ARGS} -t

norm.so: norm.cpp
	${CXX} ${CXX_FLAGS} ${SHARED_LIB_FLAGS} norm.cpp -o norm.so

clean:
	rm -rf *.pyc *.so *.png *.dot __pycache__

cleaner:
	rm -rf *.pyc *.so *.png *.dot __pycache__ resid.cpp mgU_mgR_.cpp
