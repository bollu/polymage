#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost -ansi-alias
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize
SHARED_LIB_FLAGS=-fPIC -shared
PYTHON=python3

# app inputs
CLASS='W'  # NAS-MG Problem Class
RUNS=1

APP_ARGS=--class=$(CLASS) --runs=$(RUNS)
OPT_ARGS=--pool_alloc --optimize_storage --early_free
TOOL_ARGS=--graph-gen --runs=$(RUNS) --timer
CXX_ARGS=--cxx="$(CXX)" --cxx_flags="$(CXX_FLAGS)"

all: norm.so new

new:
	$(PYTHON) main.py --mode='new' $(APP_ARGS) $(OPT_ARGS) $(TOOL_ARGS) $(CXX_ARGS)

tune:
	$(PYTHON) main.py --mode='tune' $(APP_ARGS) $(OPT_ARGS) $(TOOL_ARGS) $(CXX_ARGS)

existing:
	$(PYTHON) main.py --mode='existing' $(APP_ARGS) $(TOOL_ARGS) $(CXX_ARGS)

ready:
	$(PYTHON) main.py --mode='ready' $(APP_ARGS) $(TOOL_ARGS)

norm.so: norm.cpp
	$(CXX_ARGS) $(SHARED_LIB_FLAGS) norm.cpp -o norm.so

clean:
	rm -rf *.pyc *.so *.png *.dot __pycache__

cleaner:
	rm -rf *.pyc *.so *.png *.dot __pycache__ resid.cpp nas*.cpp
