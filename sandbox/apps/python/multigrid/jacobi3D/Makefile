
#CXX='icpc'
#CXX_FLAGS='-openmp -ipo -O3 -xhost -ansi-alias'
CXX='g++'
CXX_FLAGS='-fopenmp -O3 -march=native -ftree-vectorize'

CYCLE='V'
L=3
SIZE=31
NIT=10
NU1=10
NU2=0
NUC=0
RUNS=1

PROB_ARGS=--cycle=${CYCLE} --nit=${NIT} --runs=${RUNS} --n=${SIZE} --L=${L} --nu1=${NU1} --nuc=${NUC} --nu2=${NU2}

all: new

new:
	python3.4 main.py  --mode='new' ${PROB_ARGS} --cxx=${CXX} --cxx_flags=${CXX_FLAGS} --t

tune:
	python3.4 main.py  --mode='tune' ${PROB_ARGS} --cxx=${CXX} --cxx_flags=${CXX_FLAGS} --t

existing:
	python3.4 main.py  --mode='existing' ${PROB_ARGS} --cxx=${CXX} --cxx_flags=${CXX_FLAGS} --t

ready:
	python3.4 main.py  --mode='ready' ${PROB_ARGS}

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *cycle.cpp
