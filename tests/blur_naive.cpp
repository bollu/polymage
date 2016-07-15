#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
void pipeline_blur(int C, int R, float * input, float *& blury)
{
  #pragma omp parallel for schedule(static)
  for (int _T_i1 = 0; (_T_i1 <= (R / 16)); _T_i1 = (_T_i1 + 1))
  {
    /* users : ['blurx'] */
    float blurx[3][16][18];
    for (int _T_i2 = 0; (_T_i2 <= ((C - 1) / 16)); _T_i2 = (_T_i2 + 1))
    {
      for (int _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        int _ct0 = ((R < ((16 * _T_i1) + 15))? R: ((16 * _T_i1) + 15));
        int _ct1 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
        for (int _i1 = _ct1; (_i1 <= _ct0); _i1 = (_i1 + 1))
        {
          int _ct2 = ((C < ((16 * _T_i2) + 16))? C: ((16 * _T_i2) + 16));
          #pragma ivdep
          for (int _i2 = ((16 * _T_i2) + 1); (_i2 <= _ct2); _i2 = (_i2 + 1))
          {
            blurx[_i0][((-16 * _T_i1) + _i1)][((-16 * _T_i2) + _i2)] = ((input[(((_i0 * ((2 + R) * (2 + C))) + ((-1 + _i1) * (2 + C))) + _i2)] + input[(((_i0 * ((2 + R) * (2 + C))) + (_i1 * (2 + C))) + _i2)]) + (input[(((_i0 * ((2 + R) * (2 + C))) + ((1 + _i1) * (2 + C))) + _i2)] / 3.0));
          }
        }
      }
      for (int _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        int _ct3 = (((R - 1) < ((16 * _T_i1) + 15))? (R - 1): ((16 * _T_i1) + 15));
        int _ct4 = ((2 > (16 * _T_i1))? 2: (16 * _T_i1));
        for (int _i1 = _ct4; (_i1 <= _ct3); _i1 = (_i1 + 1))
        {
          int _ct5 = (((C - 1) < ((16 * _T_i2) + 17))? (C - 1): ((16 * _T_i2) + 17));
          int _ct6 = ((2 > (16 * _T_i2))? 2: (16 * _T_i2));
          #pragma ivdep
          for (int _i2 = _ct6; (_i2 <= _ct5); _i2 = (_i2 + 1))
          {
            blury[(((_i0 * ((-2 + R) * (-2 + C))) + ((_i1 - 2) * (-2 + C))) + (_i2 - 2))] = ((blurx[_i0][((-16 * _T_i1) + _i1)][(-1 + ((-16 * _T_i2) + _i2))] + blurx[_i0][((-16 * _T_i1) + _i1)][((-16 * _T_i2) + _i2)]) + (blurx[_i0][((-16 * _T_i1) + _i1)][(1 + ((-16 * _T_i2) + _i2))] / 3.0));
          }
        }
      }
    }
  }
}