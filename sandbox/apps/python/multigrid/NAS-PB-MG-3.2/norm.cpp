#include <stdio.h>
#include <algorithm>
#include <cmath>

using namespace std;

extern "C" void norm2u3(int N, void * rnm2_void_arg, void * rnmu_void_arg, void * r_void_arg)
{
  double * rnm2;
  rnm2 = (double *) (rnm2_void_arg);
  double * rnmu;
  rnmu = (double *) (rnmu_void_arg);
  double * r;
  r = (double *) (r_void_arg);

  long int n;
  double s, rnmu_;

  n = (N-2)*(N-2)*(N-2);
  s = 0.0;
  rnmu_ = 0.0;

  #pragma omp parallel for schedule(static) reduction(+:s) reduction(max:rnmu_)
  for(int i3 = 1; i3 < N-1; i3++){
  for(int i2 = 1; i2 < N-1; i2++){
  for(int i1 = 1; i1 < N-1; i1++){
    double a;
    double data = r[i3 * (N * N) + i2 * (N) + i1];
    s += data * data;
    a = std::fabs(data);
    rnmu_ = std::max(rnmu_, a);
  }}}

  rnm2[0] = (std::sqrt(s / n));
  rnmu[0] = rnmu_;
}
