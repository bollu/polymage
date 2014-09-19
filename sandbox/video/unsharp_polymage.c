#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>

#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
extern "C" void unsharp_mask(int  C, int  R, float thresh, float weight, const void * input_void, void * mask_void)
{
    float *input = (float *)input_void;
    float *mask = (float *)mask_void;

  //mask = (float *) (malloc((sizeof(float ) * ((3 * R) * C))));
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 <= ((R + 1) / 64)); _T_i1 = (_T_i1 + 1))
  {
    float  blurx[3][64][70];
    float  blury[3][64][70];
    float  sharpen[3][64][70];
    for (int  _T_i2 = -1; (_T_i2 <= ((C + 3) / 64)); _T_i2 = (_T_i2 + 1))
    {
      for (int  _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        int  _ct19 = (((R + 1) < ((64 * _T_i1) + 63))? (R + 1): ((64 * _T_i1) + 63));
        int  _ct20 = ((2 > (64 * _T_i1))? 2: (64 * _T_i1));
        for (int  _i1 = _ct20; (_i1 <= _ct19); _i1 = (_i1 + 1))
        {
          int  _ct21 = (((C + 3) < ((64 * _T_i2) + 69))? (C + 3): ((64 * _T_i2) + 69));
          int  _ct22 = ((0 > (64 * _T_i2))? 0: (64 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct22; (_i2 <= _ct21); _i2 = (_i2 + 1))
          {
            blurx[_i0][((-64 * _T_i1) + _i1)][((-64 * _T_i2) + _i2)] = (((((input[(((_i0 * ((R + 4) * (C + 4))) + ((-2 + _i1) * (C + 4))) + _i2)] * 0.0625f) + (input[(((_i0 * ((R + 4) * (C + 4))) + ((-1 + _i1) * (C + 4))) + _i2)] * 0.25f)) + (input[(((_i0 * ((R + 4) * (C + 4))) + (_i1 * (C + 4))) + _i2)] * 0.375f)) + (input[(((_i0 * ((R + 4) * (C + 4))) + ((1 + _i1) * (C + 4))) + _i2)] * 0.25f)) + (input[(((_i0 * ((R + 4) * (C + 4))) + ((2 + _i1) * (C + 4))) + _i2)] * 0.0625f));
          }
        }
      }
      for (int  _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        int  _ct23 = (((R + 1) < ((64 * _T_i1) + 63))? (R + 1): ((64 * _T_i1) + 63));
        int  _ct24 = ((2 > (64 * _T_i1))? 2: (64 * _T_i1));
        for (int  _i1 = _ct24; (_i1 <= _ct23); _i1 = (_i1 + 1))
        {
          int  _ct25 = (((C + 1) < ((64 * _T_i2) + 68))? (C + 1): ((64 * _T_i2) + 68));
          int  _ct26 = ((2 > ((64 * _T_i2) + 1))? 2: ((64 * _T_i2) + 1));
          #pragma ivdep
          for (int  _i2 = _ct26; (_i2 <= _ct25); _i2 = (_i2 + 1))
          {
            blury[_i0][((-64 * _T_i1) + _i1)][((-64 * _T_i2) + _i2)] = (((((blurx[_i0][((-64 * _T_i1) + _i1)][(-2 + ((-64 * _T_i2) + _i2))] * 0.0625f) + (blurx[_i0][((-64 * _T_i1) + _i1)][(-1 + ((-64 * _T_i2) + _i2))] * 0.25f)) + (blurx[_i0][((-64 * _T_i1) + _i1)][((-64 * _T_i2) + _i2)] * 0.375f)) + (blurx[_i0][((-64 * _T_i1) + _i1)][(1 + ((-64 * _T_i2) + _i2))] * 0.25f)) + (blurx[_i0][((-64 * _T_i1) + _i1)][(2 + ((-64 * _T_i2) + _i2))] * 0.0625f));
          }
        }
      }
      for (int  _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        int  _ct27 = (((R + 1) < ((64 * _T_i1) + 63))? (R + 1): ((64 * _T_i1) + 63));
        int  _ct28 = ((2 > (64 * _T_i1))? 2: (64 * _T_i1));
        for (int  _i1 = _ct28; (_i1 <= _ct27); _i1 = (_i1 + 1))
        {
          int  _ct29 = (((C + 1) < ((64 * _T_i2) + 67))? (C + 1): ((64 * _T_i2) + 67));
          int  _ct30 = ((2 > ((64 * _T_i2) + 2))? 2: ((64 * _T_i2) + 2));
          #pragma ivdep
          for (int  _i2 = _ct30; (_i2 <= _ct29); _i2 = (_i2 + 1))
          {
            sharpen[_i0][((-64 * _T_i1) + _i1)][((-64 * _T_i2) + _i2)] = ((input[(((_i0 * ((R + 4) * (C + 4))) + (_i1 * (C + 4))) + _i2)] * (1 + weight)) + (blury[_i0][((-64 * _T_i1) + _i1)][((-64 * _T_i2) + _i2)] * -(weight)));
          }
        }
      }
      for (int  _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        int  _ct31 = (((R + 1) < ((64 * _T_i1) + 63))? (R + 1): ((64 * _T_i1) + 63));
        int  _ct32 = ((2 > (64 * _T_i1))? 2: (64 * _T_i1));
        for (int  _i1 = _ct32; (_i1 <= _ct31); _i1 = (_i1 + 1))
        {
          int  _ct33 = (((C + 1) < ((64 * _T_i2) + 66))? (C + 1): ((64 * _T_i2) + 66));
          int  _ct34 = ((2 > ((64 * _T_i2) + 3))? 2: ((64 * _T_i2) + 3));
          #pragma ivdep
          for (int  _i2 = _ct34; (_i2 <= _ct33); _i2 = (_i2 + 1))
          {
            float  _ct35 = input[(((_i0 * ((R + 4) * (C + 4))) + (_i1 * (C + 4))) + _i2)];
            float  _ct36 = sharpen[_i0][((-64 * _T_i1) + _i1)][((-64 * _T_i2) + _i2)];
            float  _ct37 = ((std::abs((input[(((_i0 * ((R + 4) * (C + 4))) + (_i1 * (C + 4))) + _i2)] - blury[_i0][((-64 * _T_i1) + _i1)][((-64 * _T_i2) + _i2)])) < thresh)? _ct35: _ct36);
            mask[(((_i0 * (R * C)) + ((_i1 - 2) * C)) + (_i2 - 2))] = _ct37;
          }
        }
      }
    }
  }
}
