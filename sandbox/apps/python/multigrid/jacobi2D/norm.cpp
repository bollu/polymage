#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>

using namespace std;

extern "C" void pipeline_norm(int N, void *f_void_arg, void *uex_void_arg, void *u_void_arg, void *err_void_arg, void * rnorm_void_arg)
{
  double * f;
  f = (double *) (f_void_arg);
  double * uex;
  uex = (double *) (uex_void_arg);
  double * u;
  u = (double *) (u_void_arg);
  double * err;
  err = (double *) (err_void_arg);
  double * rnorm;
  rnorm = (double *) (rnorm_void_arg);

  double h = 1.0/(N+1);
  double invhh = 1.0/(h*h);

  double sigma_d;
  double sigma_r;

  sigma_d = 0.0;
  sigma_r = 0.0;

  #pragma omp parallel for schedule(static) reduction(+:sigma_r,sigma_d)
  for (int  i0 = 1; (i0 <= N); i0 = (i0 + 1)) {
  for (int  i1 = 1; (i1 <= N); i1 = (i1 + 1)) {
    double tmp1 = ( u[(i0-1)*(N+2) + (i1  )]
                  + u[(i0+1)*(N+2) + (i1  )]
                  + u[(i0  )*(N+2) + (i1-1)]
                  + u[(i0  )*(N+2) + (i1+1)]
                  - u[(i0  )*(N+2) + (i1  )] * 6.0
                  ) * invhh
                  - f[(i0  )*(N+2) + (i1  )];

    double tmp2 =   u[(i0  )*(N+2) + (i1  )]
                - uex[(i0  )*(N+2) + (i1  )];

    sigma_r += tmp1 * tmp1;
    sigma_d += tmp2 * tmp2;
  }}

  rnorm[0] = (std::sqrt(sigma_r) / N);
  err[0] = (std::sqrt(sigma_d) / N);
}
