CXX=icpc
CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared -DTIME
#CXX=g++
#CXX_FLAGS=-fopenmp -O3 -march=native -ftree-vectorize -fPIC -shared

all: new

new:
	python3.4 main.py  --mode='new' --cycle='V' --nit=10 --cycle_name='Vcycle' --t --runs=1 --n=31 --L=3 --nu1=10 -- nuc=0 --nu2=0

tune:
	python3.4 main.py  --mode='tune' --cycle='V' --nit=10 --cycle_name='Vcycle' --t --runs=1 --n=31 --L=3 --nu1=10 -- nuc=0 --nu2=0

existing:
	python3.4 main.py  --mode='existing' --cycle='V' --nit=10 --cycle_name='Vcycle' --t --runs=1 --n=31 --L=3 --nu1=10 -- nuc=0 --nu2=0

ready:
	python3.4 main.py  --mode='ready' --cycle='V' --nit=10 --cycle_name='Vcycle' --t --runs=1 --n=31 --L=3 --nu1=10 -- nuc=0 --nu2=0

clean:
	rm -rf *.pyc *.so *.dot *.png __pycache__

cleaner:
	rm -rf *.pyc *.so *.dot *.png __pycache__ *cycle.cpp
