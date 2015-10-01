#include <sys/time.h>

void timer__(timeval *t)
{
    gettimeofday(t, NULL);
}

void __timer(timeval *t1, timeval *t2, const char *s)
{
    gettimeofday(t2, NULL);
    printf("%s time = %g ms\n", s, (t2->tv_sec - t1->tv_sec) * 1000.0 + (t2->tv_usec - t1->tv_usec) / 1000.0);
}
