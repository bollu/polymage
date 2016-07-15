#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost
CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize
PYTHON=python3
ROOT="../../../../"

# app inputs
IMG=$(ROOT)"/../images/bayer_raw.png"
KELVIN=3700
CONTRAST=50
GAMMA=2.0
RUNS=1

APP_ARGS=--img=$(IMG) --colour_temp=$(KELVIN) --contrast=$(CONTRAST) --gamma=$(GAMMA)
OPT_ARGS=--pool_alloc --optimize_storage --early_free
TOOL_ARGS=--graph-gen --display --runs=$(RUNS) --timer
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
	rm -rf *.pyc *.so *.dot *.png __pycache__ *.cpp
