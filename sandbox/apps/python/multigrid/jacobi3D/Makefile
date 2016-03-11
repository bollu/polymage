#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost -ansi-alias
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize

PYTHON=python3.4

RUNS=1

# Multigrid specs :
CYCLE='V'  # MG Cycle : V / W
L=3  # MG levels
SIZE=31  # coarse-grid size along each dimension
NIT=10  # V / W Cycle iterations
NU1=10  # pre-smoothing steps
NU2=0  # post-smoothing steps
NUC=0  # corase-smoothing steps

PROB_ARGS=--cycle=${CYCLE} --nit=${NIT} --runs=${RUNS} --n=${SIZE} --L=${L} --nu1=${NU1} --nuc=${NUC} --nu2=${NU2}

all: new

new:
	${PYTHON} main.py  --mode='new' ${PROB_ARGS} --cxx="${CXX}" --cxx_flags="${CXX_FLAGS}" -t

tune:
	${PYTHON} main.py  --mode='tune' ${PROB_ARGS} --cxx="${CXX}" --cxx_flags="${CXX_FLAGS}" -t

existing:
	${PYTHON} main.py  --mode='existing' ${PROB_ARGS} --cxx="${CXX}" --cxx_flags="${CXX_FLAGS}" -t

ready:
	${PYTHON} main.py  --mode='ready' ${PROB_ARGS} -t

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *cycle.cpp
