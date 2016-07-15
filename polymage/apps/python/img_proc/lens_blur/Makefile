#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize 
PYTHON=python3
ROOT="../../../.."
#ROOT="."

# app inputs
IMG1="$(ROOT)/../images/left.png"
IMG2="$(ROOT)/../images/right.png"
RUNS=1

APP_ARGS=--imgl=${IMG1} --imgr=${IMG2} 
#APP_ARGS=--img=${IMG} --alpha=${ALPHA} --patch_size=${PATCH_SIZE} --search_area=${SEARCH_AREA}
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
