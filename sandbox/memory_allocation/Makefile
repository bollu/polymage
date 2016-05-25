CXX=icpc
CXX_FLAGS=-O3 -ipo -xhost -ansi-alias
#CXX=g++
#CXX_FLAGS=-O3 -march=native

all: test

test: test.out
	./test.out

test.out: simple_pool_allocator_test.cpp
	$(CXX) $(CXX_FLAGS) simple_pool_allocator.cpp simple_pool_allocator_test.cpp -o test.out

clean:
	rm -rf test.out
