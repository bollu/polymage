#include <iostream>
#include "simple_pool_allocator.h"

using namespace std;

int main() {
    pool_init();
    for (int i = 0; i < 100; i++){
        long long int buf_size = sizeof(int) * 100 + (i*i*i)%17;
        cout <<"------------------" <<endl;
        cout <<"Allocating (buffer " <<i <<") :(size = " <<buf_size <<") ..." <<endl;

        int * buf = (int*) pool_allocate(buf_size);
        cout <<"Num allocations = " <<pool_num_allocs()
             <<", Total bytes = " <<pool_total_alloc() <<endl;

        cout <<"is_alloc_in_use(buffer " <<i <<") : " <<is_alloc_in_use(buf) <<endl;

        cout <<"Deallocating (buffer " <<i <<") ... " <<endl;
        pool_deallocate(buf);

        cout <<"is_alloc_in_use (buffer " <<i <<") : " <<is_alloc_in_use(buf) <<endl;
    }
    cout <<"Destroying pool ..." <<endl;
    pool_destroy();
    return 0;
}
