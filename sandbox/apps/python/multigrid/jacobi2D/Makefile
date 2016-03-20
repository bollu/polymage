#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost -ansi-alias
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize
PYTHON=python3

RUNS=1

# app inputs
CYCLE='V'  # MG Cycle : V / W
L=3  # MG levels
SIZE=31  # coarse-grid size along each dimension
NIT=10  # V / W Cycle iterations
NU1=10  # pre-smoothing steps
NU2=0  # post-smoothing steps
NUC=0  # corase-smoothing steps

APP_ARGS=--cycle=${CYCLE} --nit=${NIT} --runs=${RUNS} --n=${SIZE} --L=${L} --nu1=${NU1} --nuc=${NUC} --nu2=${NU2}
OPT_ARGS=--pool_alloc --optimize_storage --early_free
TOOL_ARGS=--graph-gen --runs=$(RUNS) --timer
CXX_ARGS=--cxx="$(CXX)" --cxx_flags="$(CXX_FLAGS)"

all: new

new:
	$(PYTHON) main.py --mode='new' $(APP_ARGS) $(OPT_ARGS) $(TOOL_ARGS) $(CXX_ARGS)

tune:
	$(PYTHON) main.py --mode='tune' $(APP_ARGS) $(OPT_ARGS) $(TOOL_ARGS) $(CXX_ARGS)

existing:
	$(PYTHON) main.py --mode='existing' $(APP_ARGS) $(TOOL_ARGS) $(CXX_ARGS)

ready:
	$(PYTHON) main.py --mode='ready' $(APP_ARGS) $(TOOL_ARGS)

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *cycle.cpp
