#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize
PYTHON=python3
ROOT="../../../../"

# app inputs
IMG=$(ROOT)"/../images/venice_wikimedia.jpg"
RUNS=1

APP_ARGS=--img=${IMG}
OPT_ARGS=--pool_alloc --optimize_storage --early_free
TOOL_ARGS=--graph-gen --runs=${RUNS} --display --timer
CXX_ARGS=--cxx="${CXX}" --cxx_flags="${CXX_FLAGS}"

all: new

new:
	$(PYTHON) main.py --mode='new' ${CXX_ARGS} $(OPT_ARGS) ${APP_ARGS} ${TOOL_ARGS}

tune:
	$(PYTHON) main.py --mode='tune' ${CXX_ARGS} $(OPT_ARGS) ${APP_ARGS} ${TOOL_ARGS}

existing:
	$(PYTHON) main.py --mode='existing' ${CXX_ARGS} ${APP_ARGS} ${TOOL_ARGS}

ready:
	$(PYTHON) main.py --mode='ready' ${APP_ARGS} ${TOOL_ARGS}

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *.cpp
