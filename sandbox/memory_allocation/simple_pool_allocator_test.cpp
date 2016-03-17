#include<stdio.h>
#include "simple_pool_allocator.h"
int main() {
    pool_init();
    for (int i = 0; i < 100; i++){
        int * buf1 = (int*) pool_allocate(sizeof(int) * 100 + i%15);
        printf("Iteration=%d, Num allocations = %d, Total bytes = %d\n", i, pool_num_allocs(), pool_total_alloc());
        pool_deallocate(buf1);
    }
    pool_destroy();
    return 0;
}
