#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
extern "C" void  pipeline_laplacian(int  C, int  R, float  alpha, float  beta, void * img_colour_void_arg, void * laplacian_void_arg)
{
  unsigned char * img_colour_orig;
  img_colour_orig = (unsigned char *) (img_colour_void_arg);
  short unsigned int * laplacian;
  laplacian = (short unsigned int *) (laplacian_void_arg);

  float * img;
  img = (float *) (malloc((sizeof(float ) * (R * C))));
  float * remapLUT;
  remapLUT = (float *) (malloc((sizeof(float ) * 1536)));
  float * gPyramid_L0;
  gPyramid_L0 = (float *) (malloc((sizeof(float ) * ((4 * R) * C))));
  float * D_inGPyramid_L1;
  D_inGPyramid_L1 = (float *) (malloc((sizeof(float ) * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1)))));
  float * D_gPyramid_L1;
  D_gPyramid_L1 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 2) - 2) - 1) + 1)) * ((((C / 2) - 2) - 1) + 1)))));
  float * D_inGPyramid_L2;
  D_inGPyramid_L2 = (float *) (malloc((sizeof(float ) * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1)))));
  float * D_gPyramid_L2;
  D_gPyramid_L2 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 4) - 2) - 1) + 1)) * ((((C / 4) - 2) - 1) + 1)))));
  float * Dx_inGPyramid_L3;
  Dx_inGPyramid_L3 = (float *) (malloc((sizeof(float ) * (((((R / 8) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1)))));
  float * D_inGPyramid_L3;
  D_inGPyramid_L3 = (float *) (malloc((sizeof(float ) * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1)))));
  float * D_gPyramid_L3;
  D_gPyramid_L3 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 8) - 2) - 1) + 1)) * ((((C / 8) - 2) - 1) + 1)))));
  float * Dx_inGPyramid_L4;
  Dx_inGPyramid_L4 = (float *) (malloc((sizeof(float ) * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1)))));
  float * D_inGPyramid_L4;
  D_inGPyramid_L4 = (float *) (malloc((sizeof(float ) * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1)))));
  float * Dx_gPyramid_L4;
  Dx_gPyramid_L4 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 16) - 2) - 1) + 1)) * ((((C / 8) - 2) - 1) + 1)))));
  float * D_gPyramid_L4;
  D_gPyramid_L4 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 16) - 2) - 1) + 1)) * ((((C / 16) - 2) - 1) + 1)))));
  float * Ux_lPyramid_L3;
  Ux_lPyramid_L3 = (float *) (malloc((sizeof(float ) * ((4 * (((((R / 8) - 8) + 2) - 3) + 1)) * (((((C / 16) - 4) + 2) - 1) + 1)))));
  float * outLPyramid_L4;
  outLPyramid_L4 = (float *) (malloc((sizeof(float ) * ((((((R / 16) - 4) + 2) - 1) + 1) * (((((C / 16) - 4) + 2) - 1) + 1)))));
  float * Ux_outGPyramid_L3;
  Ux_outGPyramid_L3 = (float *) (malloc((sizeof(float ) * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1)))));
  float * Ux_outGPyramid_L2;
  Ux_outGPyramid_L2 = (float *) (malloc((sizeof(float ) * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1)))));
  float * outGPyramid_L2;
  outGPyramid_L2 = (float *) (malloc((sizeof(float ) * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 4) - 16) + 2) - 7) + 1)))));
  float * outGPyramid_L1;
  outGPyramid_L1 = (float *) (malloc((sizeof(float ) * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 2) - 32) + 2) - 15) + 1)))));
  float * result_ref_gray;
  result_ref_gray = (float *) (malloc((sizeof(float ) * ((-92 + R) * (-92 + C)))));

  short unsigned int *img_colour;
  img_colour = (short unsigned int *) (malloc((sizeof(short unsigned int) * (R * C * 3))));

  int off_left = 31;
  int total_pad = 92;

  int _R = R-total_pad;
  int _C = C-total_pad;

  for (int _i0 = off_left; _i0 < _R+off_left; _i0++)
  {
    for (int _i1 = off_left; _i1 < _C+off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = ((short unsigned int)(img_colour_orig[(_i0-off_left)*(_C)*3 + (_i1-off_left)*3 + _i2])) * 256;
      }
    }
  }
  for (int _i0 = off_left; _i0 < _R+off_left; _i0++)
  {
    for (int _i1 = 0; _i1 < off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_i0) * C * 3 + off_left * 3 + _i2];
      }
    }
  }
  for (int _i0 = off_left; _i0 < _R+off_left; _i0++)
  {
    for (int _i1 = _C+off_left; _i1 < C; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_i0) * C * 3 + (_C+(off_left-1)) * 3 + _i2];
      }
    }
  }

  for (int _i0 = 0; _i0 < off_left; _i0++)
  {
    for (int _i1 = 0; _i1 < off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[off_left * C * 3 + off_left * 3 + _i2];
      }
    }
  }
  for (int _i0 = 0; _i0 < off_left; _i0++)
  {
    for (int _i1 = off_left; _i1 < _C+off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[off_left * C * 3 + _i1 * 3 + _i2];
      }
    }
  }
  for (int _i0 = 0; _i0 < off_left; _i0++)
  {
    for (int _i1 = _C+off_left; _i1 < C; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[off_left * C * 3 + (_C+(off_left-1)) * 3 + _i2];
      }
    }
  }

  for (int _i0 = _R+off_left; _i0 < R; _i0++)
  {
    for (int _i1 = 0; _i1 < off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_R+(off_left-1)) * C * 3 + off_left * 3 + _i2];
      }
    }
  }
  for (int _i0 = _R+off_left; _i0 < R; _i0++)
  {
    for (int _i1 = off_left; _i1 < _C+off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_R+(off_left-1)) * C + _i1 * 3 + _i2];
      }
    }
  }
  for (int _i0 = _R+off_left; _i0 < R; _i0++)
  {
    for (int _i1 = _C+off_left; _i1 < C; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_R+(off_left-1)) * C * 3 + (_C+(off_left-1)) * 3 + _i2];
      }
    }
  }
  for (int  _i0 = -768; (_i0 <= 767); _i0 = (_i0 + 1))
  {
    remapLUT[(_i0 - -768)] = ((alpha * ((float ) (_i0) / 256.0f)) * std::exp(((-(((float ) (_i0) / 256.0f)) * ((float ) (_i0) / 256.0f)) / 2.0f)));
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 0; (_i1 < R); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 0; (_i2 < C); _i2 = (_i2 + 1))
    {
      img[((_i1 * C) + _i2)] = ((((0.299f * img_colour[(((_i1 * C * 3)) + _i2 * 3 + 2)]) + (0.587f * img_colour[(((_i1 * C * 3)) + _i2 * 3 + 1)])) + (0.114f * img_colour[(((_i1 * C * 3)) + _i2 * 3)])) / 65535.0f);
    }
  }

  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 <= ((R - 4) / 32)); _T_i1 = (_T_i1 + 1))
  {
    float  Dx_inGPyramid_L1[16][259];
    for (int  _T_i2 = -1; (_T_i2 <= ((C - 2) / 256)); _T_i2 = (_T_i2 + 1))
    {
      int  _ct20087 = ((((16 * _T_i1) + 15) < ((R / 2) - 2))? ((16 * _T_i1) + 15): ((R / 2) - 2));
      int  _ct20088 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
      for (int  _i1 = _ct20088; (_i1 <= _ct20087); _i1 = (_i1 + 1))
      {
        int  _ct20089 = (((C - 2) < ((256 * _T_i2) + 258))? (C - 2): ((256 * _T_i2) + 258));
        int  _ct20090 = ((1 > (256 * _T_i2))? 1: (256 * _T_i2));
        #pragma ivdep
        for (int  _i2 = _ct20090; (_i2 <= _ct20089); _i2 = (_i2 + 1))
        {
          Dx_inGPyramid_L1[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = ((((img[(((-1 + (2 * _i1)) * C) + _i2)] + (3.0f * img[(((2 * _i1) * C) + _i2)])) + (3.0f * img[(((1 + (2 * _i1)) * C) + _i2)])) + img[(((2 + (2 * _i1)) * C) + _i2)]) * 0.125f);
        }
      }
      if ((_T_i2 >= 0))
      {
        int  _ct20091 = ((((16 * _T_i1) + 15) < ((R / 2) - 2))? ((16 * _T_i1) + 15): ((R / 2) - 2));
        int  _ct20092 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
        for (int  _i1 = _ct20092; (_i1 <= _ct20091); _i1 = (_i1 + 1))
        {
          int  _ct20093 = (((C - 4) < ((256 * _T_i2) + 256))? (C - 4): ((256 * _T_i2) + 256));
          #pragma ivdep
          for (int  _i2 = ((256 * _T_i2) + 2); (_i2 <= _ct20093); _i2 = (_i2 + 2))
          {
            D_inGPyramid_L1[(((_i1 - 1) * ((((C / 2) - 2) - 1) + 1)) + ((_i2 / 2) - 1))] = ((((Dx_inGPyramid_L1[((-16 * _T_i1) + _i1)][(-1 + (2 * ((_i2 / 2) - (128 * _T_i2))))] + (3.0f * Dx_inGPyramid_L1[((-16 * _T_i1) + _i1)][(2 * ((_i2 / 2) - (128 * _T_i2)))])) + (3.0f * Dx_inGPyramid_L1[((-16 * _T_i1) + _i1)][(1 + (2 * ((_i2 / 2) - (128 * _T_i2))))])) + Dx_inGPyramid_L1[((-16 * _T_i1) + _i1)][(2 + (2 * ((_i2 / 2) - (128 * _T_i2))))]) * 0.125f);
          }
        }
      }
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 0; (_i1 < R); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 0; (_i2 < C); _i2 = (_i2 + 1))
      {
        int  _ct20094 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f)));
        int  _ct20095 = 0;
        int  _ct20096 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f))) > 0)? _ct20094: _ct20095);
        int  _ct20097 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f)));
        int  _ct20098 = 0;
        int  _ct20099 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f))) > 0)? _ct20097: _ct20098);
        int  _ct20100 = _ct20099;
        int  _ct20101 = 768;
        int  _ct20102 = ((_ct20096 < 768)? _ct20100: _ct20101);
        gPyramid_L0[(((_i0 * (R * C)) + (_i1 * C)) + _i2)] = (((beta * (img[((_i1 * C) + _i2)] - (_i0 * 0.333333333333f))) + (_i0 * 0.333333333333f)) + remapLUT[((_ct20102 - (256 * _i0)) - -768)]);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 <= ((R - 8) / 64)); _T_i1 = (_T_i1 + 1))
  {
    float  Dx_inGPyramid_L2[16][259];
    for (int  _T_i2 = -1; (_T_i2 <= ((C - 4) / 512)); _T_i2 = (_T_i2 + 1))
    {
      int  _ct20103 = ((((16 * _T_i1) + 15) < ((R / 4) - 2))? ((16 * _T_i1) + 15): ((R / 4) - 2));
      int  _ct20104 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
      for (int  _i1 = _ct20104; (_i1 <= _ct20103); _i1 = (_i1 + 1))
      {
        int  _ct20105 = ((((256 * _T_i2) + 258) < ((C / 2) - 2))? ((256 * _T_i2) + 258): ((C / 2) - 2));
        int  _ct20106 = ((1 > (256 * _T_i2))? 1: (256 * _T_i2));
        #pragma ivdep
        for (int  _i2 = _ct20106; (_i2 <= _ct20105); _i2 = (_i2 + 1))
        {
          Dx_inGPyramid_L2[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = ((((D_inGPyramid_L1[(((-2 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] + (3.0f * D_inGPyramid_L1[(((-1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))])) + (3.0f * D_inGPyramid_L1[(((2 * _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))])) + D_inGPyramid_L1[(((1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))]) * 0.125f);
        }
      }
      if ((_T_i2 >= 0))
      {
        int  _ct20107 = ((((16 * _T_i1) + 15) < ((R / 4) - 2))? ((16 * _T_i1) + 15): ((R / 4) - 2));
        int  _ct20108 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
        for (int  _i1 = _ct20108; (_i1 <= _ct20107); _i1 = (_i1 + 1))
        {
          int  _ct20109 = ((((256 * _T_i2) + 256) < ((C / 2) - 4))? ((256 * _T_i2) + 256): ((C / 2) - 4));
          #pragma ivdep
          for (int  _i2 = ((256 * _T_i2) + 2); (_i2 <= _ct20109); _i2 = (_i2 + 2))
          {
            D_inGPyramid_L2[(((_i1 - 1) * ((((C / 4) - 2) - 1) + 1)) + ((_i2 / 2) - 1))] = ((((Dx_inGPyramid_L2[((-16 * _T_i1) + _i1)][(-1 + (2 * ((_i2 / 2) - (128 * _T_i2))))] + (3.0f * Dx_inGPyramid_L2[((-16 * _T_i1) + _i1)][(2 * ((_i2 / 2) - (128 * _T_i2)))])) + (3.0f * Dx_inGPyramid_L2[((-16 * _T_i1) + _i1)][(1 + (2 * ((_i2 / 2) - (128 * _T_i2))))])) + Dx_inGPyramid_L2[((-16 * _T_i1) + _i1)][(2 + (2 * ((_i2 / 2) - (128 * _T_i2))))]) * 0.125f);
          }
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 <= ((R - 4) / 32)); _T_i1 = (_T_i1 + 1))
  {
    float  Dx_gPyramid_L1[4][16][259];
    for (int  _T_i2 = -1; (_T_i2 <= ((C - 2) / 256)); _T_i2 = (_T_i2 + 1))
    {
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20110 = ((((16 * _T_i1) + 15) < ((R / 2) - 2))? ((16 * _T_i1) + 15): ((R / 2) - 2));
        int  _ct20111 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
        for (int  _i1 = _ct20111; (_i1 <= _ct20110); _i1 = (_i1 + 1))
        {
          int  _ct20112 = (((C - 2) < ((256 * _T_i2) + 258))? (C - 2): ((256 * _T_i2) + 258));
          int  _ct20113 = ((1 > (256 * _T_i2))? 1: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20113; (_i2 <= _ct20112); _i2 = (_i2 + 1))
          {
            Dx_gPyramid_L1[_i0][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = ((((gPyramid_L0[(((_i0 * (R * C)) + ((-1 + (2 * _i1)) * C)) + _i2)] + (3.0f * gPyramid_L0[(((_i0 * (R * C)) + ((2 * _i1) * C)) + _i2)])) + (3.0f * gPyramid_L0[(((_i0 * (R * C)) + ((1 + (2 * _i1)) * C)) + _i2)])) + gPyramid_L0[(((_i0 * (R * C)) + ((2 + (2 * _i1)) * C)) + _i2)]) * 0.125f);
          }
        }
      }
      if ((_T_i2 >= 0))
      {
        for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
        {
          int  _ct20114 = ((((16 * _T_i1) + 15) < ((R / 2) - 2))? ((16 * _T_i1) + 15): ((R / 2) - 2));
          int  _ct20115 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
          for (int  _i1 = _ct20115; (_i1 <= _ct20114); _i1 = (_i1 + 1))
          {
            int  _ct20116 = (((C - 4) < ((256 * _T_i2) + 256))? (C - 4): ((256 * _T_i2) + 256));
            #pragma ivdep
            for (int  _i2 = ((256 * _T_i2) + 2); (_i2 <= _ct20116); _i2 = (_i2 + 2))
            {
              D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 2) - 2) - 1) + 1))) + ((_i2 / 2) - 1))] = ((((Dx_gPyramid_L1[_i0][((-16 * _T_i1) + _i1)][(-1 + (2 * ((_i2 / 2) - (128 * _T_i2))))] + (3.0f * Dx_gPyramid_L1[_i0][((-16 * _T_i1) + _i1)][(2 * ((_i2 / 2) - (128 * _T_i2)))])) + (3.0f * Dx_gPyramid_L1[_i0][((-16 * _T_i1) + _i1)][(1 + (2 * ((_i2 / 2) - (128 * _T_i2))))])) + Dx_gPyramid_L1[_i0][((-16 * _T_i1) + _i1)][(2 + (2 * ((_i2 / 2) - (128 * _T_i2))))]) * 0.125f);
            }
          }
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 <= ((R - 8) / 64)); _T_i1 = (_T_i1 + 1))
  {
    float  Dx_gPyramid_L2[4][16][259];
    for (int  _T_i2 = -1; (_T_i2 <= ((C - 4) / 512)); _T_i2 = (_T_i2 + 1))
    {
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20117 = ((((16 * _T_i1) + 15) < ((R / 4) - 2))? ((16 * _T_i1) + 15): ((R / 4) - 2));
        int  _ct20118 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
        for (int  _i1 = _ct20118; (_i1 <= _ct20117); _i1 = (_i1 + 1))
        {
          int  _ct20119 = ((((256 * _T_i2) + 258) < ((C / 2) - 2))? ((256 * _T_i2) + 258): ((C / 2) - 2));
          int  _ct20120 = ((1 > (256 * _T_i2))? 1: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20120; (_i2 <= _ct20119); _i2 = (_i2 + 1))
          {
            Dx_gPyramid_L2[_i0][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = ((((D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-2 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))] + (3.0f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))])) + (3.0f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((2 * _i1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))])) + D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))]) * 0.125f);
          }
        }
      }
      if ((_T_i2 >= 0))
      {
        for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
        {
          int  _ct20121 = ((((16 * _T_i1) + 15) < ((R / 4) - 2))? ((16 * _T_i1) + 15): ((R / 4) - 2));
          int  _ct20122 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
          for (int  _i1 = _ct20122; (_i1 <= _ct20121); _i1 = (_i1 + 1))
          {
            int  _ct20123 = ((((256 * _T_i2) + 256) < ((C / 2) - 4))? ((256 * _T_i2) + 256): ((C / 2) - 4));
            #pragma ivdep
            for (int  _i2 = ((256 * _T_i2) + 2); (_i2 <= _ct20123); _i2 = (_i2 + 2))
            {
              D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 4) - 2) - 1) + 1))) + ((_i2 / 2) - 1))] = ((((Dx_gPyramid_L2[_i0][((-16 * _T_i1) + _i1)][(-1 + (2 * ((_i2 / 2) - (128 * _T_i2))))] + (3.0f * Dx_gPyramid_L2[_i0][((-16 * _T_i1) + _i1)][(2 * ((_i2 / 2) - (128 * _T_i2)))])) + (3.0f * Dx_gPyramid_L2[_i0][((-16 * _T_i1) + _i1)][(1 + (2 * ((_i2 / 2) - (128 * _T_i2))))])) + Dx_gPyramid_L2[_i0][((-16 * _T_i1) + _i1)][(2 + (2 * ((_i2 / 2) - (128 * _T_i2))))]) * 0.125f);
            }
          }
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 1; (_i1 < ((R / 8) - 1)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 1; (_i2 < ((C / 4) - 1)); _i2 = (_i2 + 1))
    {
      Dx_inGPyramid_L3[(((_i1 - 1) * ((((C / 4) - 2) - 1) + 1)) + (_i2 - 1))] = ((((D_inGPyramid_L2[(((-2 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] + (3.0f * D_inGPyramid_L2[(((-1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))])) + (3.0f * D_inGPyramid_L2[(((2 * _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))])) + D_inGPyramid_L2[(((1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))]) * 0.125f);
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 <= ((R - 16) / 128)); _T_i1 = (_T_i1 + 1))
  {
    float  Dx_gPyramid_L3[4][16][259];
    for (int  _T_i2 = -1; (_T_i2 <= ((C - 8) / 1024)); _T_i2 = (_T_i2 + 1))
    {
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20124 = ((((16 * _T_i1) + 15) < ((R / 8) - 2))? ((16 * _T_i1) + 15): ((R / 8) - 2));
        int  _ct20125 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
        for (int  _i1 = _ct20125; (_i1 <= _ct20124); _i1 = (_i1 + 1))
        {
          int  _ct20126 = ((((256 * _T_i2) + 258) < ((C / 4) - 2))? ((256 * _T_i2) + 258): ((C / 4) - 2));
          int  _ct20127 = ((1 > (256 * _T_i2))? 1: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20127; (_i2 <= _ct20126); _i2 = (_i2 + 1))
          {
            Dx_gPyramid_L3[_i0][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = ((((D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-2 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))] + (3.0f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))])) + (3.0f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((2 * _i1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))])) + D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))]) * 0.125f);
          }
        }
      }
      if ((_T_i2 >= 0))
      {
        for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
        {
          int  _ct20128 = ((((16 * _T_i1) + 15) < ((R / 8) - 2))? ((16 * _T_i1) + 15): ((R / 8) - 2));
          int  _ct20129 = ((1 > (16 * _T_i1))? 1: (16 * _T_i1));
          for (int  _i1 = _ct20129; (_i1 <= _ct20128); _i1 = (_i1 + 1))
          {
            int  _ct20130 = ((((256 * _T_i2) + 256) < ((C / 4) - 4))? ((256 * _T_i2) + 256): ((C / 4) - 4));
            #pragma ivdep
            for (int  _i2 = ((256 * _T_i2) + 2); (_i2 <= _ct20130); _i2 = (_i2 + 2))
            {
              D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 8) - 2) - 1) + 1))) + ((_i2 / 2) - 1))] = ((((Dx_gPyramid_L3[_i0][((-16 * _T_i1) + _i1)][(-1 + (2 * ((_i2 / 2) - (128 * _T_i2))))] + (3.0f * Dx_gPyramid_L3[_i0][((-16 * _T_i1) + _i1)][(2 * ((_i2 / 2) - (128 * _T_i2)))])) + (3.0f * Dx_gPyramid_L3[_i0][((-16 * _T_i1) + _i1)][(1 + (2 * ((_i2 / 2) - (128 * _T_i2))))])) + Dx_gPyramid_L3[_i0][((-16 * _T_i1) + _i1)][(2 + (2 * ((_i2 / 2) - (128 * _T_i2))))]) * 0.125f);
            }
          }
        }
      }
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 8) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 8) - 1)); _i2 = (_i2 + 1))
    {
      D_inGPyramid_L3[(((_i1 - 1) * ((((C / 8) - 2) - 1) + 1)) + (_i2 - 1))] = ((((Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-2 + (2 * _i2)))] + (3.0f * Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + (2 * _i2)))])) + (3.0f * Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (2 * _i2))])) + Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (1 + (2 * _i2)))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 8) - 1)); _i2 = (_i2 + 1))
      {
        Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 8) - 2) - 1) + 1))) + (_i2 - 1))] = ((((D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-2 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))] + (3.0f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))])) + (3.0f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((2 * _i1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))])) + D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))]) * 0.125f);
      }
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 8) - 1)); _i2 = (_i2 + 1))
    {
      Dx_inGPyramid_L4[(((_i1 - 1) * ((((C / 8) - 2) - 1) + 1)) + (_i2 - 1))] = ((((D_inGPyramid_L3[(((-2 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] + (3.0f * D_inGPyramid_L3[(((-1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))])) + (3.0f * D_inGPyramid_L3[(((2 * _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))])) + D_inGPyramid_L3[(((1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
    {
      for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
      {
        D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 16) - 2) - 1) + 1))) + (_i2 - 1))] = ((((Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (-2 + (2 * _i2)))] + (3.0f * Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + (2 * _i2)))])) + (3.0f * Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (2 * _i2))])) + Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (1 + (2 * _i2)))]) * 0.125f);
      }
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
    {
      D_inGPyramid_L4[(((_i1 - 1) * ((((C / 16) - 2) - 1) + 1)) + (_i2 - 1))] = ((((Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-2 + (2 * _i2)))] + (3.0f * Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + (2 * _i2)))])) + (3.0f * Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (2 * _i2))])) + Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (1 + (2 * _i2)))]) * 0.125f);
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
    {
      int  _ct20131 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20132 = 0;
      int  _ct20133 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20131: _ct20132);
      int  _ct20134 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20135 = 0;
      int  _ct20136 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20134: _ct20135);
      int  _ct20137 = _ct20136;
      int  _ct20138 = 2;
      int  _ct20139 = ((_ct20133 < 2)? _ct20137: _ct20138);
      int  _ct20140 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20141 = 0;
      int  _ct20142 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20140: _ct20141);
      int  _ct20143 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20144 = 0;
      int  _ct20145 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20143: _ct20144);
      int  _ct20146 = _ct20145;
      int  _ct20147 = 2;
      int  _ct20148 = ((_ct20142 < 2)? _ct20146: _ct20147);
      int  _ct20149 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20150 = 0;
      int  _ct20151 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20149: _ct20150);
      int  _ct20152 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20153 = 0;
      int  _ct20154 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20152: _ct20153);
      int  _ct20155 = _ct20154;
      int  _ct20156 = 2;
      int  _ct20157 = ((_ct20151 < 2)? _ct20155: _ct20156);
      int  _ct20158 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20159 = 0;
      int  _ct20160 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20158: _ct20159);
      int  _ct20161 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct20162 = 0;
      int  _ct20163 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20161: _ct20162);
      int  _ct20164 = _ct20163;
      int  _ct20165 = 2;
      int  _ct20166 = ((_ct20160 < 2)? _ct20164: _ct20165);
      outLPyramid_L4[(((_i1 - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (_i2 - 1))] = (((1.0f - ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20139)) * D_gPyramid_L4[(((_ct20148 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]) + (((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20157) * D_gPyramid_L4[((((_ct20166 + 1) * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]));
    }
  }
  for (int  _i1 = 3; (_i1 < ((R / 8) - 5)); _i1 = (_i1 + 2))
  {
    for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
    {
      Ux_outGPyramid_L3[(((_i1 - 3) * (((((C / 16) - 4) + 2) - 1) + 1)) + (_i2 - 1))] = ((0.25f * outLPyramid_L4[(((((_i1 / 2) + 1) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]) + (0.75f * outLPyramid_L4[((((_i1 / 2) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]));
    }
    if ((R >= ((8 * _i1) + 56)))
    {
      for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
      {
        Ux_outGPyramid_L3[((((_i1 + 1) - 3) * (((((C / 16) - 4) + 2) - 1) + 1)) + (_i2 - 1))] = ((0.25f * outLPyramid_L4[((((((_i1 + 1) / 2) - 1) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]) + (0.75f * outLPyramid_L4[(((((_i1 + 1) / 2) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]));
      }
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 3; (_i1 < ((R / 8) - 5)); _i1 = (_i1 + 2))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
      {
        Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((_i1 - 3) * (((((C / 16) - 4) + 2) - 1) + 1))) + (_i2 - 1))] = ((0.25f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]));
      }
      if ((R >= ((8 * _i1) + 56)))
      {
        #pragma ivdep
        for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
        {
          Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + (((_i1 + 1) - 3) * (((((C / 16) - 4) + 2) - 1) + 1))) + (_i2 - 1))] = ((0.25f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + (((((_i1 + 1) / 2) - 1) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((((_i1 + 1) / 2) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]));
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = -1; (_T_i1 < ((R + 16) / 64)); _T_i1 = (_T_i1 + 1))
  {
    float  U_lPyramid_L3[4][14][256];
    float  outLPyramid_L3[14][256];
    float  outGPyramid_L3[14][256];
    for (int  _T_i2 = 0; (_T_i2 <= ((C - 48) / 2048)); _T_i2 = (_T_i2 + 1))
    {
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20167 = ((((16 * _T_i1) + 26) < ((R / 4) - 12))? ((16 * _T_i1) + 26): ((R / 4) - 12));
        int  _ct20168 = ((6 > (16 * _T_i1))? 6: (16 * _T_i1));
        for (int  _i1 = _ct20168; (_i1 <= _ct20167); _i1 = (_i1 + 2))
        {
          int  _ct20169 = ((((256 * _T_i2) + 255) < ((C / 8) - 6))? ((256 * _T_i2) + 255): ((C / 8) - 6));
          int  _ct20170 = ((3 > (256 * _T_i2))? 3: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20170; (_i2 <= _ct20169); _i2 = (_i2 + 1))
          {
            float  _ct20171 = ((0.25f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1))) + (((_i2 / 2) - 1) - 1))]) + (0.75f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((_i2 / 2) - 1))]));
            float  _ct20172 = ((0.25f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1))) + (((_i2 / 2) + 1) - 1))]) + (0.75f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((_i2 / 2) - 1))]));
            float  _ct20173 = (((_i2 % 2) == 0)? _ct20171: _ct20172);
            U_lPyramid_L3[_i0][((_i1 / 2) - (8 * _T_i1))][((-256 * _T_i2) + _i2)] = (D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))] - _ct20173);
          }
        }
      }
      int  _ct20174 = ((((16 * _T_i1) + 24) < ((R / 4) - 12))? ((16 * _T_i1) + 24): ((R / 4) - 12));
      int  _ct20175 = ((6 > ((16 * _T_i1) + 2))? 6: ((16 * _T_i1) + 2));
      for (int  _i1 = _ct20175; (_i1 <= _ct20174); _i1 = (_i1 + 2))
      {
        int  _ct20176 = ((((256 * _T_i2) + 255) < ((C / 8) - 6))? ((256 * _T_i2) + 255): ((C / 8) - 6));
        int  _ct20177 = ((3 > (256 * _T_i2))? 3: (256 * _T_i2));
        #pragma ivdep
        for (int  _i2 = _ct20177; (_i2 <= _ct20176); _i2 = (_i2 + 1))
        {
          int  _ct20178 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20179 = 0;
          int  _ct20180 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20178: _ct20179);
          int  _ct20181 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20182 = 0;
          int  _ct20183 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20181: _ct20182);
          int  _ct20184 = _ct20183;
          int  _ct20185 = 2;
          int  _ct20186 = ((_ct20180 < 2)? _ct20184: _ct20185);
          int  _ct20187 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20188 = 0;
          int  _ct20189 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20187: _ct20188);
          int  _ct20190 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20191 = 0;
          int  _ct20192 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20190: _ct20191);
          int  _ct20193 = _ct20192;
          int  _ct20194 = 2;
          int  _ct20195 = ((_ct20189 < 2)? _ct20193: _ct20194);
          int  _ct20196 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20197 = 0;
          int  _ct20198 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20196: _ct20197);
          int  _ct20199 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20200 = 0;
          int  _ct20201 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20199: _ct20200);
          int  _ct20202 = _ct20201;
          int  _ct20203 = 2;
          int  _ct20204 = ((_ct20198 < 2)? _ct20202: _ct20203);
          int  _ct20205 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20206 = 0;
          int  _ct20207 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20205: _ct20206);
          int  _ct20208 = (int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20209 = 0;
          int  _ct20210 = (((int ) ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20208: _ct20209);
          int  _ct20211 = _ct20210;
          int  _ct20212 = 2;
          int  _ct20213 = ((_ct20207 < 2)? _ct20211: _ct20212);
          outLPyramid_L3[((_i1 / 2) - (8 * _T_i1))][((-256 * _T_i2) + _i2)] = (((1.0f - ((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20186)) * U_lPyramid_L3[_ct20195][((_i1 / 2) - (8 * _T_i1))][((-256 * _T_i2) + _i2)]) + (((D_inGPyramid_L3[(((-1 + (_i1 / 2)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20204) * U_lPyramid_L3[(_ct20213 + 1)][((_i1 / 2) - (8 * _T_i1))][((-256 * _T_i2) + _i2)]));
        }
      }
      int  _ct20214 = ((((16 * _T_i1) + 22) < ((R / 4) - 12))? ((16 * _T_i1) + 22): ((R / 4) - 12));
      int  _ct20215 = ((6 > ((16 * _T_i1) + 4))? 6: ((16 * _T_i1) + 4));
      for (int  _i1 = _ct20215; (_i1 <= _ct20214); _i1 = (_i1 + 2))
      {
        int  _ct20216 = ((((256 * _T_i2) + 255) < ((C / 8) - 6))? ((256 * _T_i2) + 255): ((C / 8) - 6));
        int  _ct20217 = ((3 > (256 * _T_i2))? 3: (256 * _T_i2));
        #pragma ivdep
        for (int  _i2 = _ct20217; (_i2 <= _ct20216); _i2 = (_i2 + 1))
        {
          float  _ct20218 = ((0.25f * Ux_outGPyramid_L3[(((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1)) + (((_i2 / 2) - 1) - 1))]) + (0.75f * Ux_outGPyramid_L3[(((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1)) + ((_i2 / 2) - 1))]));
          float  _ct20219 = ((0.25f * Ux_outGPyramid_L3[(((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1)) + (((_i2 / 2) + 1) - 1))]) + (0.75f * Ux_outGPyramid_L3[(((-3 + (_i1 / 2)) * (((((C / 16) - 4) + 2) - 1) + 1)) + ((_i2 / 2) - 1))]));
          float  _ct20220 = (((_i2 % 2) == 0)? _ct20218: _ct20219);
          outGPyramid_L3[((_i1 / 2) - (8 * _T_i1))][((-256 * _T_i2) + _i2)] = (outLPyramid_L3[((_i1 / 2) - (8 * _T_i1))][((-256 * _T_i2) + _i2)] + _ct20220);
        }
      }
      int  _ct20221 = ((((16 * _T_i1) + 20) < ((R / 4) - 14))? ((16 * _T_i1) + 20): ((R / 4) - 14));
      int  _ct20222 = ((8 > ((16 * _T_i1) + 6))? 8: ((16 * _T_i1) + 6));
      for (int  _i1 = _ct20222; (_i1 <= _ct20221); _i1 = (_i1 + 2))
      {
        int  _ct20223 = ((((256 * _T_i2) + 255) < ((C / 8) - 6))? ((256 * _T_i2) + 255): ((C / 8) - 6));
        int  _ct20224 = ((3 > (256 * _T_i2))? 3: (256 * _T_i2));
        #pragma ivdep
        for (int  _i2 = _ct20224; (_i2 <= _ct20223); _i2 = (_i2 + 1))
        {
          Ux_outGPyramid_L2[(((_i1 - 7) * (((((C / 8) - 8) + 2) - 3) + 1)) + (_i2 - 3))] = ((0.25f * outGPyramid_L3[((((-16 * _T_i1) + _i1) / 2) - 1)][((-256 * _T_i2) + _i2)]) + (0.75f * outGPyramid_L3[(((-16 * _T_i1) + _i1) / 2)][((-256 * _T_i2) + _i2)]));
        }
      }
      if ((_T_i1 >= 0))
      {
        int  _ct20225 = ((((16 * _T_i1) + 21) < ((R / 4) - 14))? ((16 * _T_i1) + 21): ((R / 4) - 14));
        for (int  _i1 = ((16 * _T_i1) + 7); (_i1 <= _ct20225); _i1 = (_i1 + 2))
        {
          int  _ct20226 = ((((256 * _T_i2) + 255) < ((C / 8) - 6))? ((256 * _T_i2) + 255): ((C / 8) - 6));
          int  _ct20227 = ((3 > (256 * _T_i2))? 3: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20227; (_i2 <= _ct20226); _i2 = (_i2 + 1))
          {
            Ux_outGPyramid_L2[(((_i1 - 7) * (((((C / 8) - 8) + 2) - 3) + 1)) + (_i2 - 3))] = ((0.25f * outGPyramid_L3[((((-16 * _T_i1) + _i1) / 2) + 1)][((-256 * _T_i2) + _i2)]) + (0.75f * outGPyramid_L3[(((-16 * _T_i1) + _i1) / 2)][((-256 * _T_i2) + _i2)]));
          }
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 < ((R + 8) / 64)); _T_i1 = (_T_i1 + 1))
  {
    float  Ux_lPyramid_L2[4][16][131];
    float  U_lPyramid_L2[4][16][262];
    float  outLPyramid_L2[16][262];
    for (int  _T_i2 = 0; (_T_i2 <= ((C - 48) / 1024)); _T_i2 = (_T_i2 + 1))
    {
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20228 = ((((16 * _T_i1) + 14) < ((R / 4) - 14))? ((16 * _T_i1) + 14): ((R / 4) - 14));
        int  _ct20229 = ((8 > (16 * _T_i1))? 8: (16 * _T_i1));
        for (int  _i1 = _ct20229; (_i1 <= _ct20228); _i1 = (_i1 + 2))
        {
          int  _ct20230 = ((((256 * _T_i2) + 260) < ((C / 4) - 12))? ((256 * _T_i2) + 260): ((C / 4) - 12));
          int  _ct20231 = ((6 > (256 * _T_i2))? 6: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20231; (_i2 <= _ct20230); _i2 = (_i2 + 2))
          {
            Ux_lPyramid_L2[_i0][((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((((_i1 / 2) - 1) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]) + (0.75f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]));
          }
        }
      }
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20232 = ((((16 * _T_i1) + 15) < ((R / 4) - 14))? ((16 * _T_i1) + 15): ((R / 4) - 14));
        int  _ct20233 = ((7 > ((16 * _T_i1) + 1))? 7: ((16 * _T_i1) + 1));
        for (int  _i1 = _ct20233; (_i1 <= _ct20232); _i1 = (_i1 + 2))
        {
          int  _ct20234 = ((((256 * _T_i2) + 260) < ((C / 4) - 12))? ((256 * _T_i2) + 260): ((C / 4) - 12));
          int  _ct20235 = ((6 > (256 * _T_i2))? 6: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20235; (_i2 <= _ct20234); _i2 = (_i2 + 2))
          {
            Ux_lPyramid_L2[_i0][((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]) + (0.75f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]));
          }
        }
      }
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20236 = ((((16 * _T_i1) + 15) < ((R / 4) - 14))? ((16 * _T_i1) + 15): ((R / 4) - 14));
        int  _ct20237 = ((7 > (16 * _T_i1))? 7: (16 * _T_i1));
        for (int  _i1 = _ct20237; (_i1 <= _ct20236); _i1 = (_i1 + 1))
        {
          int  _ct20238 = ((((256 * _T_i2) + 260) < ((C / 4) - 14))? ((256 * _T_i2) + 260): ((C / 4) - 14));
          int  _ct20239 = ((7 > ((256 * _T_i2) + 1))? 7: ((256 * _T_i2) + 1));
          #pragma ivdep
          for (int  _i2 = _ct20239; (_i2 <= _ct20238); _i2 = (_i2 + 1))
          {
            float  _ct20240 = ((0.25f * Ux_lPyramid_L2[_i0][((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) - 1)]) + (0.75f * Ux_lPyramid_L2[_i0][((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
            float  _ct20241 = ((0.25f * Ux_lPyramid_L2[_i0][((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) + 1)]) + (0.75f * Ux_lPyramid_L2[_i0][((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
            float  _ct20242 = (((_i2 % 2) == 0)? _ct20240: _ct20241);
            U_lPyramid_L2[_i0][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = (D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))] - _ct20242);
          }
        }
      }
      int  _ct20243 = ((((16 * _T_i1) + 15) < ((R / 4) - 14))? ((16 * _T_i1) + 15): ((R / 4) - 14));
      int  _ct20244 = ((7 > (16 * _T_i1))? 7: (16 * _T_i1));
      for (int  _i1 = _ct20244; (_i1 <= _ct20243); _i1 = (_i1 + 1))
      {
        int  _ct20245 = ((((256 * _T_i2) + 259) < ((C / 4) - 14))? ((256 * _T_i2) + 259): ((C / 4) - 14));
        int  _ct20246 = ((7 > ((256 * _T_i2) + 2))? 7: ((256 * _T_i2) + 2));
        #pragma ivdep
        for (int  _i2 = _ct20246; (_i2 <= _ct20245); _i2 = (_i2 + 1))
        {
          int  _ct20247 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20248 = 0;
          int  _ct20249 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20247: _ct20248);
          int  _ct20250 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20251 = 0;
          int  _ct20252 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20250: _ct20251);
          int  _ct20253 = _ct20252;
          int  _ct20254 = 2;
          int  _ct20255 = ((_ct20249 < 2)? _ct20253: _ct20254);
          int  _ct20256 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20257 = 0;
          int  _ct20258 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20256: _ct20257);
          int  _ct20259 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20260 = 0;
          int  _ct20261 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20259: _ct20260);
          int  _ct20262 = _ct20261;
          int  _ct20263 = 2;
          int  _ct20264 = ((_ct20258 < 2)? _ct20262: _ct20263);
          int  _ct20265 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20266 = 0;
          int  _ct20267 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20265: _ct20266);
          int  _ct20268 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20269 = 0;
          int  _ct20270 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20268: _ct20269);
          int  _ct20271 = _ct20270;
          int  _ct20272 = 2;
          int  _ct20273 = ((_ct20267 < 2)? _ct20271: _ct20272);
          int  _ct20274 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20275 = 0;
          int  _ct20276 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20274: _ct20275);
          int  _ct20277 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20278 = 0;
          int  _ct20279 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20277: _ct20278);
          int  _ct20280 = _ct20279;
          int  _ct20281 = 2;
          int  _ct20282 = ((_ct20276 < 2)? _ct20280: _ct20281);
          outLPyramid_L2[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = (((1.0f - ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20255)) * U_lPyramid_L2[_ct20264][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)]) + (((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20273) * U_lPyramid_L2[(_ct20282 + 1)][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)]));
        }
      }
      int  _ct20283 = ((((16 * _T_i1) + 15) < ((R / 4) - 14))? ((16 * _T_i1) + 15): ((R / 4) - 14));
      int  _ct20284 = ((7 > (16 * _T_i1))? 7: (16 * _T_i1));
      for (int  _i1 = _ct20284; (_i1 <= _ct20283); _i1 = (_i1 + 1))
      {
        int  _ct20285 = ((((256 * _T_i2) + 258) < ((C / 4) - 14))? ((256 * _T_i2) + 258): ((C / 4) - 14));
        int  _ct20286 = ((7 > ((256 * _T_i2) + 3))? 7: ((256 * _T_i2) + 3));
        #pragma ivdep
        for (int  _i2 = _ct20286; (_i2 <= _ct20285); _i2 = (_i2 + 1))
        {
          float  _ct20287 = ((0.25f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + (((_i2 / 2) - 1) - 3))]) + (0.75f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + ((_i2 / 2) - 3))]));
          float  _ct20288 = ((0.25f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + (((_i2 / 2) + 1) - 3))]) + (0.75f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + ((_i2 / 2) - 3))]));
          float  _ct20289 = (((_i2 % 2) == 0)? _ct20287: _ct20288);
          outGPyramid_L2[(((_i1 - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (_i2 - 7))] = (outLPyramid_L2[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] + _ct20289);
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 0; (_T_i1 < (((R + 4) / 32) - 1)); _T_i1 = (_T_i1 + 1))
  {
    float  Ux_lPyramid_L1[4][16][131];
    float  Ux_outGPyramid_L1[16][131];
    float  U_lPyramid_L1[4][16][262];
    float  outLPyramid_L1[16][262];
    for (int  _T_i2 = 0; (_T_i2 <= ((C - 56) / 512)); _T_i2 = (_T_i2 + 1))
    {
      if ((_T_i1 >= 1))
      {
        for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
        {
          int  _ct20290 = ((((16 * _T_i1) + 14) < ((R / 2) - 30))? ((16 * _T_i1) + 14): ((R / 2) - 30));
          for (int  _i1 = (16 * _T_i1); (_i1 <= _ct20290); _i1 = (_i1 + 2))
          {
            int  _ct20291 = ((((256 * _T_i2) + 260) < ((C / 2) - 28))? ((256 * _T_i2) + 260): ((C / 2) - 28));
            int  _ct20292 = ((14 > (256 * _T_i2))? 14: (256 * _T_i2));
            #pragma ivdep
            for (int  _i2 = _ct20292; (_i2 <= _ct20291); _i2 = (_i2 + 2))
            {
              Ux_lPyramid_L1[_i0][((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((((_i1 / 2) - 1) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]) + (0.75f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]));
            }
          }
        }
      }
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20293 = ((((16 * _T_i1) + 15) < ((R / 2) - 30))? ((16 * _T_i1) + 15): ((R / 2) - 30));
        int  _ct20294 = ((15 > ((16 * _T_i1) + 1))? 15: ((16 * _T_i1) + 1));
        for (int  _i1 = _ct20294; (_i1 <= _ct20293); _i1 = (_i1 + 2))
        {
          int  _ct20295 = ((((256 * _T_i2) + 260) < ((C / 2) - 28))? ((256 * _T_i2) + 260): ((C / 2) - 28));
          int  _ct20296 = ((14 > (256 * _T_i2))? 14: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20296; (_i2 <= _ct20295); _i2 = (_i2 + 2))
          {
            Ux_lPyramid_L1[_i0][((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]) + (0.75f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]));
          }
        }
      }
      if ((_T_i1 >= 1))
      {
        int  _ct20297 = ((((16 * _T_i1) + 14) < ((R / 2) - 30))? ((16 * _T_i1) + 14): ((R / 2) - 30));
        for (int  _i1 = (16 * _T_i1); (_i1 <= _ct20297); _i1 = (_i1 + 2))
        {
          int  _ct20298 = ((((256 * _T_i2) + 260) < ((C / 2) - 28))? ((256 * _T_i2) + 260): ((C / 2) - 28));
          int  _ct20299 = ((14 > (256 * _T_i2))? 14: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20299; (_i2 <= _ct20298); _i2 = (_i2 + 2))
          {
            Ux_outGPyramid_L1[((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * outGPyramid_L2[(((((_i1 / 2) - 1) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + (_i2 / 2)))]) + (0.75f * outGPyramid_L2[((((_i1 / 2) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + (_i2 / 2)))]));
          }
        }
      }
      int  _ct20300 = ((((16 * _T_i1) + 15) < ((R / 2) - 30))? ((16 * _T_i1) + 15): ((R / 2) - 30));
      int  _ct20301 = ((15 > ((16 * _T_i1) + 1))? 15: ((16 * _T_i1) + 1));
      for (int  _i1 = _ct20301; (_i1 <= _ct20300); _i1 = (_i1 + 2))
      {
        int  _ct20302 = ((((256 * _T_i2) + 260) < ((C / 2) - 28))? ((256 * _T_i2) + 260): ((C / 2) - 28));
        int  _ct20303 = ((14 > (256 * _T_i2))? 14: (256 * _T_i2));
        #pragma ivdep
        for (int  _i2 = _ct20303; (_i2 <= _ct20302); _i2 = (_i2 + 2))
        {
          Ux_outGPyramid_L1[((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * outGPyramid_L2[(((((_i1 / 2) + 1) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + (_i2 / 2)))]) + (0.75f * outGPyramid_L2[((((_i1 / 2) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + (_i2 / 2)))]));
        }
      }
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20304 = ((((16 * _T_i1) + 15) < ((R / 2) - 30))? ((16 * _T_i1) + 15): ((R / 2) - 30));
        int  _ct20305 = ((15 > (16 * _T_i1))? 15: (16 * _T_i1));
        for (int  _i1 = _ct20305; (_i1 <= _ct20304); _i1 = (_i1 + 1))
        {
          int  _ct20306 = ((((256 * _T_i2) + 260) < ((C / 2) - 30))? ((256 * _T_i2) + 260): ((C / 2) - 30));
          int  _ct20307 = ((15 > ((256 * _T_i2) + 1))? 15: ((256 * _T_i2) + 1));
          #pragma ivdep
          for (int  _i2 = _ct20307; (_i2 <= _ct20306); _i2 = (_i2 + 1))
          {
            float  _ct20308 = ((0.25f * Ux_lPyramid_L1[_i0][((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) - 1)]) + (0.75f * Ux_lPyramid_L1[_i0][((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
            float  _ct20309 = ((0.25f * Ux_lPyramid_L1[_i0][((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) + 1)]) + (0.75f * Ux_lPyramid_L1[_i0][((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
            float  _ct20310 = (((_i2 % 2) == 0)? _ct20308: _ct20309);
            U_lPyramid_L1[_i0][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = (D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))] - _ct20310);
          }
        }
      }
      int  _ct20311 = ((((16 * _T_i1) + 15) < ((R / 2) - 30))? ((16 * _T_i1) + 15): ((R / 2) - 30));
      int  _ct20312 = ((15 > (16 * _T_i1))? 15: (16 * _T_i1));
      for (int  _i1 = _ct20312; (_i1 <= _ct20311); _i1 = (_i1 + 1))
      {
        int  _ct20313 = ((((256 * _T_i2) + 259) < ((C / 2) - 30))? ((256 * _T_i2) + 259): ((C / 2) - 30));
        int  _ct20314 = ((15 > ((256 * _T_i2) + 2))? 15: ((256 * _T_i2) + 2));
        #pragma ivdep
        for (int  _i2 = _ct20314; (_i2 <= _ct20313); _i2 = (_i2 + 1))
        {
          int  _ct20315 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20316 = 0;
          int  _ct20317 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20315: _ct20316);
          int  _ct20318 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20319 = 0;
          int  _ct20320 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20318: _ct20319);
          int  _ct20321 = _ct20320;
          int  _ct20322 = 2;
          int  _ct20323 = ((_ct20317 < 2)? _ct20321: _ct20322);
          int  _ct20324 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20325 = 0;
          int  _ct20326 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20324: _ct20325);
          int  _ct20327 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20328 = 0;
          int  _ct20329 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20327: _ct20328);
          int  _ct20330 = _ct20329;
          int  _ct20331 = 2;
          int  _ct20332 = ((_ct20326 < 2)? _ct20330: _ct20331);
          int  _ct20333 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20334 = 0;
          int  _ct20335 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20333: _ct20334);
          int  _ct20336 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20337 = 0;
          int  _ct20338 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20336: _ct20337);
          int  _ct20339 = _ct20338;
          int  _ct20340 = 2;
          int  _ct20341 = ((_ct20335 < 2)? _ct20339: _ct20340);
          int  _ct20342 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20343 = 0;
          int  _ct20344 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20342: _ct20343);
          int  _ct20345 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
          int  _ct20346 = 0;
          int  _ct20347 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct20345: _ct20346);
          int  _ct20348 = _ct20347;
          int  _ct20349 = 2;
          int  _ct20350 = ((_ct20344 < 2)? _ct20348: _ct20349);
          outLPyramid_L1[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = (((1.0f - ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20323)) * U_lPyramid_L1[_ct20332][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)]) + (((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20341) * U_lPyramid_L1[(_ct20350 + 1)][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)]));
        }
      }
      int  _ct20351 = ((((16 * _T_i1) + 15) < ((R / 2) - 30))? ((16 * _T_i1) + 15): ((R / 2) - 30));
      int  _ct20352 = ((15 > (16 * _T_i1))? 15: (16 * _T_i1));
      for (int  _i1 = _ct20352; (_i1 <= _ct20351); _i1 = (_i1 + 1))
      {
        int  _ct20353 = ((((256 * _T_i2) + 258) < ((C / 2) - 30))? ((256 * _T_i2) + 258): ((C / 2) - 30));
        int  _ct20354 = ((15 > ((256 * _T_i2) + 3))? 15: ((256 * _T_i2) + 3));
        #pragma ivdep
        for (int  _i2 = _ct20354; (_i2 <= _ct20353); _i2 = (_i2 + 1))
        {
          float  _ct20355 = ((0.25f * Ux_outGPyramid_L1[((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) - 1)]) + (0.75f * Ux_outGPyramid_L1[((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
          float  _ct20356 = ((0.25f * Ux_outGPyramid_L1[((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) + 1)]) + (0.75f * Ux_outGPyramid_L1[((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
          float  _ct20357 = (((_i2 % 2) == 0)? _ct20355: _ct20356);
          outGPyramid_L1[(((_i1 - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (_i2 - 15))] = (outLPyramid_L1[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] + _ct20357);
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _T_i1 = 1; (_T_i1 < (((R + 2) / 16) - 3)); _T_i1 = (_T_i1 + 1))
  {
    float  Ux_lPyramid_L0[4][16][131];
    float  Ux_result_ref_gray[16][131];
    float  U_lPyramid_L0[4][16][262];
    float  outLPyramid_L0[16][262];
    for (int  _T_i2 = 0; (_T_i2 <= ((C - 60) / 256)); _T_i2 = (_T_i2 + 1))
    {
      if ((_T_i1 >= 2))
      {
        for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
        {
          int  _ct20358 = (((R - 62) < ((16 * _T_i1) + 14))? (R - 62): ((16 * _T_i1) + 14));
          for (int  _i1 = (16 * _T_i1); (_i1 <= _ct20358); _i1 = (_i1 + 2))
          {
            int  _ct20359 = (((C - 60) < ((256 * _T_i2) + 260))? (C - 60): ((256 * _T_i2) + 260));
            int  _ct20360 = ((30 > (256 * _T_i2))? 30: (256 * _T_i2));
            #pragma ivdep
            for (int  _i2 = _ct20360; (_i2 <= _ct20359); _i2 = (_i2 + 2))
            {
              Ux_lPyramid_L0[_i0][((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((((_i1 / 2) - 1) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]) + (0.75f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]));
            }
          }
        }
      }
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20361 = (((R - 62) < ((16 * _T_i1) + 15))? (R - 62): ((16 * _T_i1) + 15));
        int  _ct20362 = ((31 > ((16 * _T_i1) + 1))? 31: ((16 * _T_i1) + 1));
        for (int  _i1 = _ct20362; (_i1 <= _ct20361); _i1 = (_i1 + 2))
        {
          int  _ct20363 = (((C - 60) < ((256 * _T_i2) + 260))? (C - 60): ((256 * _T_i2) + 260));
          int  _ct20364 = ((30 > (256 * _T_i2))? 30: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20364; (_i2 <= _ct20363); _i2 = (_i2 + 2))
          {
            Ux_lPyramid_L0[_i0][((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]) + (0.75f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + (_i2 / 2)))]));
          }
        }
      }
      if ((_T_i1 >= 2))
      {
        int  _ct20365 = (((R - 62) < ((16 * _T_i1) + 14))? (R - 62): ((16 * _T_i1) + 14));
        for (int  _i1 = (16 * _T_i1); (_i1 <= _ct20365); _i1 = (_i1 + 2))
        {
          int  _ct20366 = (((C - 60) < ((256 * _T_i2) + 260))? (C - 60): ((256 * _T_i2) + 260));
          int  _ct20367 = ((30 > (256 * _T_i2))? 30: (256 * _T_i2));
          #pragma ivdep
          for (int  _i2 = _ct20367; (_i2 <= _ct20366); _i2 = (_i2 + 2))
          {
            Ux_result_ref_gray[((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * outGPyramid_L1[(((((_i1 / 2) - 1) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + (_i2 / 2)))]) + (0.75f * outGPyramid_L1[((((_i1 / 2) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + (_i2 / 2)))]));
          }
        }
      }
      int  _ct20368 = (((R - 62) < ((16 * _T_i1) + 15))? (R - 62): ((16 * _T_i1) + 15));
      int  _ct20369 = ((31 > ((16 * _T_i1) + 1))? 31: ((16 * _T_i1) + 1));
      for (int  _i1 = _ct20369; (_i1 <= _ct20368); _i1 = (_i1 + 2))
      {
        int  _ct20370 = (((C - 60) < ((256 * _T_i2) + 260))? (C - 60): ((256 * _T_i2) + 260));
        int  _ct20371 = ((30 > (256 * _T_i2))? 30: (256 * _T_i2));
        #pragma ivdep
        for (int  _i2 = _ct20371; (_i2 <= _ct20370); _i2 = (_i2 + 2))
        {
          Ux_result_ref_gray[((-16 * _T_i1) + _i1)][((_i2 / 2) - (128 * _T_i2))] = ((0.25f * outGPyramid_L1[(((((_i1 / 2) + 1) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + (_i2 / 2)))]) + (0.75f * outGPyramid_L1[((((_i1 / 2) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + (_i2 / 2)))]));
        }
      }
      for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
      {
        int  _ct20372 = (((R - 62) < ((16 * _T_i1) + 15))? (R - 62): ((16 * _T_i1) + 15));
        int  _ct20373 = ((31 > (16 * _T_i1))? 31: (16 * _T_i1));
        for (int  _i1 = _ct20373; (_i1 <= _ct20372); _i1 = (_i1 + 1))
        {
          int  _ct20374 = (((C - 62) < ((256 * _T_i2) + 260))? (C - 62): ((256 * _T_i2) + 260));
          int  _ct20375 = ((31 > ((256 * _T_i2) + 1))? 31: ((256 * _T_i2) + 1));
          #pragma ivdep
          for (int  _i2 = _ct20375; (_i2 <= _ct20374); _i2 = (_i2 + 1))
          {
            float  _ct20376 = ((0.25f * Ux_lPyramid_L0[_i0][((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) - 1)]) + (0.75f * Ux_lPyramid_L0[_i0][((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
            float  _ct20377 = ((0.25f * Ux_lPyramid_L0[_i0][((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) + 1)]) + (0.75f * Ux_lPyramid_L0[_i0][((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
            float  _ct20378 = (((_i2 % 2) == 0)? _ct20376: _ct20377);
            U_lPyramid_L0[_i0][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = (gPyramid_L0[(((_i0 * (R * C)) + (_i1 * C)) + _i2)] - _ct20378);
          }
        }
      }
      int  _ct20379 = (((R - 62) < ((16 * _T_i1) + 15))? (R - 62): ((16 * _T_i1) + 15));
      int  _ct20380 = ((31 > (16 * _T_i1))? 31: (16 * _T_i1));
      for (int  _i1 = _ct20380; (_i1 <= _ct20379); _i1 = (_i1 + 1))
      {
        int  _ct20381 = (((C - 62) < ((256 * _T_i2) + 259))? (C - 62): ((256 * _T_i2) + 259));
        int  _ct20382 = ((31 > ((256 * _T_i2) + 2))? 31: ((256 * _T_i2) + 2));
        #pragma ivdep
        for (int  _i2 = _ct20382; (_i2 <= _ct20381); _i2 = (_i2 + 1))
        {
          int  _ct20383 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20384 = 0;
          int  _ct20385 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20383: _ct20384);
          int  _ct20386 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20387 = 0;
          int  _ct20388 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20386: _ct20387);
          int  _ct20389 = _ct20388;
          int  _ct20390 = 2;
          int  _ct20391 = ((_ct20385 < 2)? _ct20389: _ct20390);
          int  _ct20392 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20393 = 0;
          int  _ct20394 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20392: _ct20393);
          int  _ct20395 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20396 = 0;
          int  _ct20397 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20395: _ct20396);
          int  _ct20398 = _ct20397;
          int  _ct20399 = 2;
          int  _ct20400 = ((_ct20394 < 2)? _ct20398: _ct20399);
          int  _ct20401 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20402 = 0;
          int  _ct20403 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20401: _ct20402);
          int  _ct20404 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20405 = 0;
          int  _ct20406 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20404: _ct20405);
          int  _ct20407 = _ct20406;
          int  _ct20408 = 2;
          int  _ct20409 = ((_ct20403 < 2)? _ct20407: _ct20408);
          int  _ct20410 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20411 = 0;
          int  _ct20412 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20410: _ct20411);
          int  _ct20413 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
          int  _ct20414 = 0;
          int  _ct20415 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct20413: _ct20414);
          int  _ct20416 = _ct20415;
          int  _ct20417 = 2;
          int  _ct20418 = ((_ct20412 < 2)? _ct20416: _ct20417);
          outLPyramid_L0[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] = (((1.0f - ((img[((_i1 * C) + _i2)] * (float ) (3)) - _ct20391)) * U_lPyramid_L0[_ct20400][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)]) + (((img[((_i1 * C) + _i2)] * (float ) (3)) - _ct20409) * U_lPyramid_L0[(_ct20418 + 1)][((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)]));
        }
      }
      int  _ct20419 = (((R - 62) < ((16 * _T_i1) + 15))? (R - 62): ((16 * _T_i1) + 15));
      int  _ct20420 = ((31 > (16 * _T_i1))? 31: (16 * _T_i1));
      for (int  _i1 = _ct20420; (_i1 <= _ct20419); _i1 = (_i1 + 1))
      {
        int  _ct20421 = (((C - 62) < ((256 * _T_i2) + 258))? (C - 62): ((256 * _T_i2) + 258));
        int  _ct20422 = ((31 > ((256 * _T_i2) + 3))? 31: ((256 * _T_i2) + 3));
        #pragma ivdep
        for (int  _i2 = _ct20422; (_i2 <= _ct20421); _i2 = (_i2 + 1))
        {
          float  _ct20423 = ((0.25f * Ux_result_ref_gray[((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) - 1)]) + (0.75f * Ux_result_ref_gray[((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
          float  _ct20424 = ((0.25f * Ux_result_ref_gray[((-16 * _T_i1) + _i1)][((((-256 * _T_i2) + _i2) / 2) + 1)]) + (0.75f * Ux_result_ref_gray[((-16 * _T_i1) + _i1)][(((-256 * _T_i2) + _i2) / 2)]));
          float  _ct20425 = (((_i2 % 2) == 0)? _ct20423: _ct20424);
          result_ref_gray[(((_i1 - 31) * (-92 + C)) + (_i2 - 31))] = (outLPyramid_L0[((-16 * _T_i1) + _i1)][((-256 * _T_i2) + _i2)] + _ct20425);
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 31; (_i1 < (R - 61)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 31; (_i2 < (C - 61)); _i2 = (_i2 + 1))
    {
      #pragma ivdep
      for (int  _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        float  _ct20426 = ((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f));
        float  _ct20427 = 0.0f;
        float  _ct20428 = ((((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f)) > 0.0f)? _ct20426: _ct20427);
        float  _ct20429 = ((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f));
        float  _ct20430 = 0.0f;
        float  _ct20431 = ((((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f)) > 0.0f)? _ct20429: _ct20430);
        float  _ct20432 = _ct20431;
        float  _ct20433 = 1.0f;
        float  _ct20434 = ((_ct20428 < 1.0f)? _ct20432: _ct20433);
        laplacian[((((_i1 - 31) * (-92 + C) * 3)) + (_i2 - 31) * 3 + _i0)] = (short unsigned int ) ((_ct20434 * 65535.0f));
      }
    }
  }
  free(img);
  free(img_colour);
  free(remapLUT);
  free(gPyramid_L0);
  free(D_inGPyramid_L1);
  free(D_gPyramid_L1);
  free(D_inGPyramid_L2);
  free(D_gPyramid_L2);
  free(Dx_inGPyramid_L3);
  free(D_inGPyramid_L3);
  free(D_gPyramid_L3);
  free(Dx_inGPyramid_L4);
  free(D_inGPyramid_L4);
  free(Dx_gPyramid_L4);
  free(D_gPyramid_L4);
  free(Ux_lPyramid_L3);
  free(outLPyramid_L4);
  free(Ux_outGPyramid_L3);
  free(Ux_outGPyramid_L2);
  free(outGPyramid_L2);
  free(outGPyramid_L1);
  free(result_ref_gray);
}

extern "C" void  pipeline_laplacian_naive(int  C, int  R, float  alpha, float  beta, void * img_colour_void_arg, void * laplacian_void_arg)
{
  unsigned char * img_colour_orig;
  img_colour_orig = (unsigned char *) (img_colour_void_arg);
  short unsigned int * laplacian;
  laplacian = (short unsigned int *) (laplacian_void_arg);

  float * img;
  img = (float *) (malloc((sizeof(float ) * (R * C))));
  float * remapLUT;
  remapLUT = (float *) (malloc((sizeof(float ) * 1536)));
  float * Dx_inGPyramid_L1;
  Dx_inGPyramid_L1 = (float *) (malloc((sizeof(float ) * (((((R / 2) - 2) - 1) + 1) * (-2 + C)))));
  float * gPyramid_L0;
  gPyramid_L0 = (float *) (malloc((sizeof(float ) * ((4 * R) * C))));
  float * D_inGPyramid_L1;
  D_inGPyramid_L1 = (float *) (malloc((sizeof(float ) * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1)))));
  float * Dx_gPyramid_L1;
  Dx_gPyramid_L1 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 2) - 2) - 1) + 1)) * (-2 + C)))));
  float * D_gPyramid_L1;
  D_gPyramid_L1 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 2) - 2) - 1) + 1)) * ((((C / 2) - 2) - 1) + 1)))));
  float * Dx_inGPyramid_L2;
  Dx_inGPyramid_L2 = (float *) (malloc((sizeof(float ) * (((((R / 4) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1)))));
  float * D_inGPyramid_L2;
  D_inGPyramid_L2 = (float *) (malloc((sizeof(float ) * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1)))));
  float * Dx_gPyramid_L2;
  Dx_gPyramid_L2 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 4) - 2) - 1) + 1)) * ((((C / 2) - 2) - 1) + 1)))));
  float * Ux_lPyramid_L0;
  Ux_lPyramid_L0 = (float *) (malloc((sizeof(float ) * ((4 * (-92 + R)) * (((((C / 2) - 32) + 2) - 15) + 1)))));
  float * D_gPyramid_L2;
  D_gPyramid_L2 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 4) - 2) - 1) + 1)) * ((((C / 4) - 2) - 1) + 1)))));
  float * Dx_inGPyramid_L3;
  Dx_inGPyramid_L3 = (float *) (malloc((sizeof(float ) * (((((R / 8) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1)))));
  float * U_lPyramid_L0;
  U_lPyramid_L0 = (float *) (malloc((sizeof(float ) * ((4 * (-92 + R)) * (-92 + C)))));
  float * D_inGPyramid_L3;
  D_inGPyramid_L3 = (float *) (malloc((sizeof(float ) * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1)))));
  float * Dx_gPyramid_L3;
  Dx_gPyramid_L3 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 8) - 2) - 1) + 1)) * ((((C / 4) - 2) - 1) + 1)))));
  float * Ux_lPyramid_L1;
  Ux_lPyramid_L1 = (float *) (malloc((sizeof(float ) * ((4 * (((((R / 2) - 32) + 2) - 15) + 1)) * (((((C / 4) - 16) + 2) - 7) + 1)))));
  float * outLPyramid_L0;
  outLPyramid_L0 = (float *) (malloc((sizeof(float ) * ((-92 + R) * (-92 + C)))));
  float * D_gPyramid_L3;
  D_gPyramid_L3 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 8) - 2) - 1) + 1)) * ((((C / 8) - 2) - 1) + 1)))));
  float * Dx_inGPyramid_L4;
  Dx_inGPyramid_L4 = (float *) (malloc((sizeof(float ) * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1)))));
  float * U_lPyramid_L1;
  U_lPyramid_L1 = (float *) (malloc((sizeof(float ) * ((4 * (((((R / 2) - 32) + 2) - 15) + 1)) * (((((C / 2) - 32) + 2) - 15) + 1)))));
  float * D_inGPyramid_L4;
  D_inGPyramid_L4 = (float *) (malloc((sizeof(float ) * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1)))));
  float * Dx_gPyramid_L4;
  Dx_gPyramid_L4 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 16) - 2) - 1) + 1)) * ((((C / 8) - 2) - 1) + 1)))));
  float * Ux_lPyramid_L2;
  Ux_lPyramid_L2 = (float *) (malloc((sizeof(float ) * ((4 * (((((R / 4) - 16) + 2) - 7) + 1)) * (((((C / 8) - 8) + 2) - 3) + 1)))));
  float * outLPyramid_L1;
  outLPyramid_L1 = (float *) (malloc((sizeof(float ) * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 2) - 32) + 2) - 15) + 1)))));
  float * D_gPyramid_L4;
  D_gPyramid_L4 = (float *) (malloc((sizeof(float ) * ((4 * ((((R / 16) - 2) - 1) + 1)) * ((((C / 16) - 2) - 1) + 1)))));
  float * U_lPyramid_L2;
  U_lPyramid_L2 = (float *) (malloc((sizeof(float ) * ((4 * (((((R / 4) - 16) + 2) - 7) + 1)) * (((((C / 4) - 16) + 2) - 7) + 1)))));
  float * Ux_lPyramid_L3;
  Ux_lPyramid_L3 = (float *) (malloc((sizeof(float ) * ((4 * (((((R / 8) - 8) + 2) - 3) + 1)) * (((((C / 16) - 4) + 2) - 1) + 1)))));
  float * outLPyramid_L2;
  outLPyramid_L2 = (float *) (malloc((sizeof(float ) * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 4) - 16) + 2) - 7) + 1)))));
  float * outLPyramid_L4;
  outLPyramid_L4 = (float *) (malloc((sizeof(float ) * ((((((R / 16) - 4) + 2) - 1) + 1) * (((((C / 16) - 4) + 2) - 1) + 1)))));
  float * U_lPyramid_L3;
  U_lPyramid_L3 = (float *) (malloc((sizeof(float ) * ((4 * (((((R / 8) - 8) + 2) - 3) + 1)) * (((((C / 8) - 8) + 2) - 3) + 1)))));
  float * Ux_outGPyramid_L3;
  Ux_outGPyramid_L3 = (float *) (malloc((sizeof(float ) * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1)))));
  float * outLPyramid_L3;
  outLPyramid_L3 = (float *) (malloc((sizeof(float ) * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 8) - 8) + 2) - 3) + 1)))));
  float * outGPyramid_L3;
  outGPyramid_L3 = (float *) (malloc((sizeof(float ) * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 8) - 8) + 2) - 3) + 1)))));
  float * Ux_outGPyramid_L2;
  Ux_outGPyramid_L2 = (float *) (malloc((sizeof(float ) * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1)))));
  float * outGPyramid_L2;
  outGPyramid_L2 = (float *) (malloc((sizeof(float ) * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 4) - 16) + 2) - 7) + 1)))));
  float * Ux_outGPyramid_L1;
  Ux_outGPyramid_L1 = (float *) (malloc((sizeof(float ) * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 4) - 16) + 2) - 7) + 1)))));
  float * outGPyramid_L1;
  outGPyramid_L1 = (float *) (malloc((sizeof(float ) * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 2) - 32) + 2) - 15) + 1)))));
  float * Ux_result_ref_gray;
  Ux_result_ref_gray = (float *) (malloc((sizeof(float ) * ((-92 + R) * (((((C / 2) - 32) + 2) - 15) + 1)))));
  float * result_ref_gray;
  result_ref_gray = (float *) (malloc((sizeof(float ) * ((-92 + R) * (-92 + C)))));

  short unsigned int *img_colour;
  img_colour = (short unsigned int *) (malloc((sizeof(short unsigned int) * (R * C * 3))));

  int off_left = 31;
  int total_pad = 92;

  int _R = R-total_pad;
  int _C = C-total_pad;

  for (int _i0 = off_left; _i0 < _R+off_left; _i0++)
  {
    for (int _i1 = off_left; _i1 < _C+off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = ((unsigned short int)(img_colour_orig[(_i0-off_left)*(_C)*3 + (_i1-off_left)*3 + _i2])) * 256;
      }
    }
  }
  for (int _i0 = off_left; _i0 < _R+off_left; _i0++)
  {
    for (int _i1 = 0; _i1 < off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_i0) * C * 3 + off_left * 3 + _i2];
      }
    }
  }
  for (int _i0 = off_left; _i0 < _R+off_left; _i0++)
  {
    for (int _i1 = _C+off_left; _i1 < C; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_i0) * C * 3 + (_C+(off_left-1)) * 3 + _i2];
      }
    }
  }

  for (int _i0 = 0; _i0 < off_left; _i0++)
  {
    for (int _i1 = 0; _i1 < off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[off_left * C * 3 + off_left * 3 + _i2];
      }
    }
  }
  for (int _i0 = 0; _i0 < off_left; _i0++)
  {
    for (int _i1 = off_left; _i1 < _C+off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[off_left * C * 3 + _i1 * 3 + _i2];
      }
    }
  }
  for (int _i0 = 0; _i0 < off_left; _i0++)
  {
    for (int _i1 = _C+off_left; _i1 < C; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[off_left * C * 3 + (_C+(off_left-1)) * 3 + _i2];
      }
    }
  }

  for (int _i0 = _R+off_left; _i0 < R; _i0++)
  {
    for (int _i1 = 0; _i1 < off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_R+(off_left-1)) * C * 3 + off_left * 3 + _i2];
      }
    }
  }
  for (int _i0 = _R+off_left; _i0 < R; _i0++)
  {
    for (int _i1 = off_left; _i1 < _C+off_left; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_R+(off_left-1)) * C + _i1 * 3 + _i2];
      }
    }
  }
  for (int _i0 = _R+off_left; _i0 < R; _i0++)
  {
    for (int _i1 = _C+off_left; _i1 < C; _i1++)
    {
      for (int _i2 = 0; _i2 <= 2; _i2++)
      {
        img_colour[_i0 * C * 3 + _i1 * 3 + _i2] = img_colour[(_R+(off_left-1)) * C * 3 + (_C+(off_left-1)) * 3 + _i2];
      }
    }
  }

  for (int  _i0 = -768; (_i0 <= 767); _i0 = (_i0 + 1))
  {
    remapLUT[(_i0 - -768)] = ((alpha * ((float ) (_i0) / 256.0f)) * std::exp(((-(((float ) (_i0) / 256.0f)) * ((float ) (_i0) / 256.0f)) / 2.0f)));
  }

  #pragma omp parallel for schedule(static)
  for (int  _i1 = 0; (_i1 < R); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 0; (_i2 < C); _i2 = (_i2 + 1))
    {
      img[((_i1 * C) + _i2)] = ((((0.299f * img_colour[(((_i1 * C * 3)) + _i2 * 3 + 2)]) + (0.587f * img_colour[(((_i1 * C * 3)) + _i2 * 3 + 1)])) + (0.114f * img_colour[(((_i1 * C * 3)) + _i2 * 3)])) / 65535.0f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 0; (_i1 < R); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 0; (_i2 < C); _i2 = (_i2 + 1))
      {
        int  _ct0 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f)));
        int  _ct1 = 0;
        int  _ct2 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f))) > 0)? _ct0: _ct1);
        int  _ct3 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f)));
        int  _ct4 = 0;
        int  _ct5 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (768.0f))) > 0)? _ct3: _ct4);
        int  _ct6 = _ct5;
        int  _ct7 = 768;
        int  _ct8 = ((_ct2 < 768)? _ct6: _ct7);
        gPyramid_L0[(((_i0 * (R * C)) + (_i1 * C)) + _i2)] = (((beta * (img[((_i1 * C) + _i2)] - (_i0 * 0.333333333333f))) + (_i0 * 0.333333333333f)) + remapLUT[((_ct8 - (256 * _i0)) - -768)]);
      }
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 2) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < (C - 1)); _i2 = (_i2 + 1))
      {
        Dx_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * (-2 + C))) + ((_i1 - 1) * (-2 + C))) + (_i2 - 1))] = ((((gPyramid_L0[(((_i0 * (R * C)) + ((-1 + (2 * _i1)) * C)) + _i2)] + (3.0f * gPyramid_L0[(((_i0 * (R * C)) + ((2 * _i1) * C)) + _i2)])) + (3.0f * gPyramid_L0[(((_i0 * (R * C)) + ((1 + (2 * _i1)) * C)) + _i2)])) + gPyramid_L0[(((_i0 * (R * C)) + ((2 + (2 * _i1)) * C)) + _i2)]) * 0.125f);
      }
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 2) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 2) - 1)); _i2 = (_i2 + 1))
      {
        D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 2) - 2) - 1) + 1))) + (_i2 - 1))] = ((((Dx_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * (-2 + C))) + ((-1 + _i1) * (-2 + C))) + (-2 + (2 * _i2)))] + (3.0f * Dx_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * (-2 + C))) + ((-1 + _i1) * (-2 + C))) + (-1 + (2 * _i2)))])) + (3.0f * Dx_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * (-2 + C))) + ((-1 + _i1) * (-2 + C))) + (2 * _i2))])) + Dx_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * (-2 + C))) + ((-1 + _i1) * (-2 + C))) + (1 + (2 * _i2)))]) * 0.125f);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 1; (_i1 < ((R / 2) - 1)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 1; (_i2 < (C - 1)); _i2 = (_i2 + 1))
    {
      Dx_inGPyramid_L1[(((_i1 - 1) * (-2 + C)) + (_i2 - 1))] = ((((img[(((-1 + (2 * _i1)) * C) + _i2)] + (3.0f * img[(((2 * _i1) * C) + _i2)])) + (3.0f * img[(((1 + (2 * _i1)) * C) + _i2)])) + img[(((2 + (2 * _i1)) * C) + _i2)]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 4) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 2) - 1)); _i2 = (_i2 + 1))
      {
        Dx_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 2) - 2) - 1) + 1))) + (_i2 - 1))] = ((((D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-2 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))] + (3.0f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))])) + (3.0f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((2 * _i1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))])) + D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))]) * 0.125f);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 1; (_i1 < ((R / 2) - 1)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 1; (_i2 < ((C / 2) - 1)); _i2 = (_i2 + 1))
    {
      D_inGPyramid_L1[(((_i1 - 1) * ((((C / 2) - 2) - 1) + 1)) + (_i2 - 1))] = ((((Dx_inGPyramid_L1[(((-1 + _i1) * (-2 + C)) + (-2 + (2 * _i2)))] + (3.0f * Dx_inGPyramid_L1[(((-1 + _i1) * (-2 + C)) + (-1 + (2 * _i2)))])) + (3.0f * Dx_inGPyramid_L1[(((-1 + _i1) * (-2 + C)) + (2 * _i2))])) + Dx_inGPyramid_L1[(((-1 + _i1) * (-2 + C)) + (1 + (2 * _i2)))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 4) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 4) - 1)); _i2 = (_i2 + 1))
      {
        D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 4) - 2) - 1) + 1))) + (_i2 - 1))] = ((((Dx_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 2) - 2) - 1) + 1))) + (-2 + (2 * _i2)))] + (3.0f * Dx_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + (2 * _i2)))])) + (3.0f * Dx_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 2) - 2) - 1) + 1))) + (2 * _i2))])) + Dx_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 2) - 2) - 1) + 1))) + (1 + (2 * _i2)))]) * 0.125f);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 1; (_i1 < ((R / 4) - 1)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 1; (_i2 < ((C / 2) - 1)); _i2 = (_i2 + 1))
    {
      Dx_inGPyramid_L2[(((_i1 - 1) * ((((C / 2) - 2) - 1) + 1)) + (_i2 - 1))] = ((((D_inGPyramid_L1[(((-2 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] + (3.0f * D_inGPyramid_L1[(((-1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))])) + (3.0f * D_inGPyramid_L1[(((2 * _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))])) + D_inGPyramid_L1[(((1 + (2 * _i1)) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 8) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 4) - 1)); _i2 = (_i2 + 1))
      {
        Dx_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 4) - 2) - 1) + 1))) + (_i2 - 1))] = ((((D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-2 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))] + (3.0f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))])) + (3.0f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((2 * _i1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))])) + D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))]) * 0.125f);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 1; (_i1 < ((R / 4) - 1)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 1; (_i2 < ((C / 4) - 1)); _i2 = (_i2 + 1))
    {
      D_inGPyramid_L2[(((_i1 - 1) * ((((C / 4) - 2) - 1) + 1)) + (_i2 - 1))] = ((((Dx_inGPyramid_L2[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-2 + (2 * _i2)))] + (3.0f * Dx_inGPyramid_L2[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + (2 * _i2)))])) + (3.0f * Dx_inGPyramid_L2[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (2 * _i2))])) + Dx_inGPyramid_L2[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (1 + (2 * _i2)))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 8) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 8) - 1)); _i2 = (_i2 + 1))
      {
        D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 8) - 2) - 1) + 1))) + (_i2 - 1))] = ((((Dx_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 4) - 2) - 1) + 1))) + (-2 + (2 * _i2)))] + (3.0f * Dx_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + (2 * _i2)))])) + (3.0f * Dx_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 4) - 2) - 1) + 1))) + (2 * _i2))])) + Dx_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 4) - 2) - 1) + 1))) + (1 + (2 * _i2)))]) * 0.125f);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 1; (_i1 < ((R / 8) - 1)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 1; (_i2 < ((C / 4) - 1)); _i2 = (_i2 + 1))
    {
      Dx_inGPyramid_L3[(((_i1 - 1) * ((((C / 4) - 2) - 1) + 1)) + (_i2 - 1))] = ((((D_inGPyramid_L2[(((-2 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] + (3.0f * D_inGPyramid_L2[(((-1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))])) + (3.0f * D_inGPyramid_L2[(((2 * _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))])) + D_inGPyramid_L2[(((1 + (2 * _i1)) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 8) - 1)); _i2 = (_i2 + 1))
      {
        Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 8) - 2) - 1) + 1))) + (_i2 - 1))] = ((((D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-2 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))] + (3.0f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))])) + (3.0f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((2 * _i1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))])) + D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))]) * 0.125f);
      }
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 8) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 8) - 1)); _i2 = (_i2 + 1))
    {
      D_inGPyramid_L3[(((_i1 - 1) * ((((C / 8) - 2) - 1) + 1)) + (_i2 - 1))] = ((((Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-2 + (2 * _i2)))] + (3.0f * Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + (2 * _i2)))])) + (3.0f * Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (2 * _i2))])) + Dx_inGPyramid_L3[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (1 + (2 * _i2)))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
    {
      for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
      {
        D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((_i1 - 1) * ((((C / 16) - 2) - 1) + 1))) + (_i2 - 1))] = ((((Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (-2 + (2 * _i2)))] + (3.0f * Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + (2 * _i2)))])) + (3.0f * Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (2 * _i2))])) + Dx_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (1 + (2 * _i2)))]) * 0.125f);
      }
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 8) - 1)); _i2 = (_i2 + 1))
    {
      Dx_inGPyramid_L4[(((_i1 - 1) * ((((C / 8) - 2) - 1) + 1)) + (_i2 - 1))] = ((((D_inGPyramid_L3[(((-2 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] + (3.0f * D_inGPyramid_L3[(((-1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))])) + (3.0f * D_inGPyramid_L3[(((2 * _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))])) + D_inGPyramid_L3[(((1 + (2 * _i1)) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 3; (_i1 < ((R / 8) - 5)); _i1 = (_i1 + 2))
    {
      #pragma ivdep
      for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
      {
        Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((_i1 - 3) * (((((C / 16) - 4) + 2) - 1) + 1))) + (_i2 - 1))] = ((0.25f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]));
      }
      if ((R >= ((8 * _i1) + 56)))
      {
        #pragma ivdep
        for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
        {
          Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + (((_i1 + 1) - 3) * (((((C / 16) - 4) + 2) - 1) + 1))) + (_i2 - 1))] = ((0.25f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + (((((_i1 + 1) / 2) - 1) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L4[(((_i0 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((((_i1 + 1) / 2) - 1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]));
        }
      }
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
    {
      D_inGPyramid_L4[(((_i1 - 1) * ((((C / 16) - 2) - 1) + 1)) + (_i2 - 1))] = ((((Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-2 + (2 * _i2)))] + (3.0f * Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + (2 * _i2)))])) + (3.0f * Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (2 * _i2))])) + Dx_inGPyramid_L4[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (1 + (2 * _i2)))]) * 0.125f);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 3; (_i1 < ((R / 8) - 5)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 3; (_i2 < ((C / 8) - 5)); _i2 = (_i2 + 1))
      {
        float  _ct9 = ((0.25f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1))) + (((_i2 / 2) - 1) - 1))]) + (0.75f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((_i2 / 2) - 1))]));
        float  _ct10 = ((0.25f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1))) + (((_i2 / 2) + 1) - 1))]) + (0.75f * Ux_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1))) + ((_i2 / 2) - 1))]));
        float  _ct11 = (((_i2 % 2) == 0)? _ct9: _ct10);
        U_lPyramid_L3[(((_i0 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((_i1 - 3) * (((((C / 8) - 8) + 2) - 3) + 1))) + (_i2 - 3))] = (D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))] - _ct11);
      }
    }
  }
  for (int  _i1 = 1; (_i1 < ((R / 16) - 1)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
    {
      int  _ct12 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct13 = 0;
      int  _ct14 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct12: _ct13);
      int  _ct15 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct16 = 0;
      int  _ct17 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct15: _ct16);
      int  _ct18 = _ct17;
      int  _ct19 = 2;
      int  _ct20 = ((_ct14 < 2)? _ct18: _ct19);
      int  _ct21 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct22 = 0;
      int  _ct23 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct21: _ct22);
      int  _ct24 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct25 = 0;
      int  _ct26 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct24: _ct25);
      int  _ct27 = _ct26;
      int  _ct28 = 2;
      int  _ct29 = ((_ct23 < 2)? _ct27: _ct28);
      int  _ct30 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct31 = 0;
      int  _ct32 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct30: _ct31);
      int  _ct33 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct34 = 0;
      int  _ct35 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct33: _ct34);
      int  _ct36 = _ct35;
      int  _ct37 = 2;
      int  _ct38 = ((_ct32 < 2)? _ct36: _ct37);
      int  _ct39 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct40 = 0;
      int  _ct41 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct39: _ct40);
      int  _ct42 = (int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct43 = 0;
      int  _ct44 = (((int ) ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct42: _ct43);
      int  _ct45 = _ct44;
      int  _ct46 = 2;
      int  _ct47 = ((_ct41 < 2)? _ct45: _ct46);
      outLPyramid_L4[(((_i1 - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (_i2 - 1))] = (((1.0f - ((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct20)) * D_gPyramid_L4[(((_ct29 * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]) + (((D_inGPyramid_L4[(((-1 + _i1) * ((((C / 16) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct38) * D_gPyramid_L4[((((_ct47 + 1) * (((((R / 16) - 2) - 1) + 1) * ((((C / 16) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 16) - 2) - 1) + 1))) + (-1 + _i2))]));
    }
  }
  for (int  _i1 = 3; (_i1 < ((R / 8) - 5)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 3; (_i2 < ((C / 8) - 5)); _i2 = (_i2 + 1))
    {
      int  _ct48 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct49 = 0;
      int  _ct50 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct48: _ct49);
      int  _ct51 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct52 = 0;
      int  _ct53 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct51: _ct52);
      int  _ct54 = _ct53;
      int  _ct55 = 2;
      int  _ct56 = ((_ct50 < 2)? _ct54: _ct55);
      int  _ct57 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct58 = 0;
      int  _ct59 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct57: _ct58);
      int  _ct60 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct61 = 0;
      int  _ct62 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct60: _ct61);
      int  _ct63 = _ct62;
      int  _ct64 = 2;
      int  _ct65 = ((_ct59 < 2)? _ct63: _ct64);
      int  _ct66 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct67 = 0;
      int  _ct68 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct66: _ct67);
      int  _ct69 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct70 = 0;
      int  _ct71 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct69: _ct70);
      int  _ct72 = _ct71;
      int  _ct73 = 2;
      int  _ct74 = ((_ct68 < 2)? _ct72: _ct73);
      int  _ct75 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct76 = 0;
      int  _ct77 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct75: _ct76);
      int  _ct78 = (int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct79 = 0;
      int  _ct80 = (((int ) ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct78: _ct79);
      int  _ct81 = _ct80;
      int  _ct82 = 2;
      int  _ct83 = ((_ct77 < 2)? _ct81: _ct82);
      outLPyramid_L3[(((_i1 - 3) * (((((C / 8) - 8) + 2) - 3) + 1)) + (_i2 - 3))] = (((1.0f - ((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct56)) * U_lPyramid_L3[(((_ct65 * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((-3 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1))) + (-3 + _i2))]) + (((D_inGPyramid_L3[(((-1 + _i1) * ((((C / 8) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct74) * U_lPyramid_L3[((((_ct83 + 1) * ((((((R / 8) - 8) + 2) - 3) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((-3 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1))) + (-3 + _i2))]));
    }
  }
  for (int  _i1 = 3; (_i1 < ((R / 8) - 5)); _i1 = (_i1 + 2))
  {
    for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
    {
      Ux_outGPyramid_L3[(((_i1 - 3) * (((((C / 16) - 4) + 2) - 1) + 1)) + (_i2 - 1))] = ((0.25f * outLPyramid_L4[(((((_i1 / 2) + 1) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]) + (0.75f * outLPyramid_L4[((((_i1 / 2) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]));
    }
    if ((R >= ((8 * _i1) + 56)))
    {
      for (int  _i2 = 1; (_i2 < ((C / 16) - 1)); _i2 = (_i2 + 1))
      {
        Ux_outGPyramid_L3[((((_i1 + 1) - 3) * (((((C / 16) - 4) + 2) - 1) + 1)) + (_i2 - 1))] = ((0.25f * outLPyramid_L4[((((((_i1 + 1) / 2) - 1) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]) + (0.75f * outLPyramid_L4[(((((_i1 + 1) / 2) - 1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (-1 + _i2))]));
      }
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 7; (_i1 < ((R / 4) - 13)); _i1 = (_i1 + 2))
    {
      #pragma ivdep
      for (int  _i2 = 3; (_i2 < ((C / 8) - 5)); _i2 = (_i2 + 1))
      {
        Ux_lPyramid_L2[(((_i0 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((_i1 - 7) * (((((C / 8) - 8) + 2) - 3) + 1))) + (_i2 - 3))] = ((0.25f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))]));
      }
      if ((R >= ((4 * _i1) + 60)))
      {
        #pragma ivdep
        for (int  _i2 = 3; (_i2 < ((C / 8) - 5)); _i2 = (_i2 + 1))
        {
          Ux_lPyramid_L2[(((_i0 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + (((_i1 + 1) - 7) * (((((C / 8) - 8) + 2) - 3) + 1))) + (_i2 - 3))] = ((0.25f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + (((((_i1 + 1) / 2) - 1) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L3[(((_i0 * (((((R / 8) - 2) - 1) + 1) * ((((C / 8) - 2) - 1) + 1))) + ((((_i1 + 1) / 2) - 1) * ((((C / 8) - 2) - 1) + 1))) + (-1 + _i2))]));
        }
      }
    }
  }
  for (int  _i1 = 3; (_i1 < ((R / 8) - 5)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 3; (_i2 < ((C / 8) - 5)); _i2 = (_i2 + 1))
    {
      float  _ct84 = ((0.25f * Ux_outGPyramid_L3[(((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (((_i2 / 2) - 1) - 1))]) + (0.75f * Ux_outGPyramid_L3[(((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1)) + ((_i2 / 2) - 1))]));
      float  _ct85 = ((0.25f * Ux_outGPyramid_L3[(((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1)) + (((_i2 / 2) + 1) - 1))]) + (0.75f * Ux_outGPyramid_L3[(((-3 + _i1) * (((((C / 16) - 4) + 2) - 1) + 1)) + ((_i2 / 2) - 1))]));
      float  _ct86 = (((_i2 % 2) == 0)? _ct84: _ct85);
      outGPyramid_L3[(((_i1 - 3) * (((((C / 8) - 8) + 2) - 3) + 1)) + (_i2 - 3))] = (outLPyramid_L3[(((-3 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + (-3 + _i2))] + _ct86);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 7; (_i1 < ((R / 4) - 13)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 7; (_i2 < ((C / 4) - 13)); _i2 = (_i2 + 1))
      {
        float  _ct87 = ((0.25f * Ux_lPyramid_L2[(((_i0 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1))) + (((_i2 / 2) - 1) - 3))]) + (0.75f * Ux_lPyramid_L2[(((_i0 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((_i2 / 2) - 3))]));
        float  _ct88 = ((0.25f * Ux_lPyramid_L2[(((_i0 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1))) + (((_i2 / 2) + 1) - 3))]) + (0.75f * Ux_lPyramid_L2[(((_i0 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1))) + ((_i2 / 2) - 3))]));
        float  _ct89 = (((_i2 % 2) == 0)? _ct87: _ct88);
        U_lPyramid_L2[(((_i0 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((_i1 - 7) * (((((C / 4) - 16) + 2) - 7) + 1))) + (_i2 - 7))] = (D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))] - _ct89);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 7; (_i1 < ((R / 4) - 13)); _i1 = (_i1 + 2))
  {
    #pragma ivdep
    for (int  _i2 = 3; (_i2 < ((C / 8) - 5)); _i2 = (_i2 + 1))
    {
      Ux_outGPyramid_L2[(((_i1 - 7) * (((((C / 8) - 8) + 2) - 3) + 1)) + (_i2 - 3))] = ((0.25f * outGPyramid_L3[(((((_i1 / 2) + 1) - 3) * (((((C / 8) - 8) + 2) - 3) + 1)) + (-3 + _i2))]) + (0.75f * outGPyramid_L3[((((_i1 / 2) - 3) * (((((C / 8) - 8) + 2) - 3) + 1)) + (-3 + _i2))]));
    }
    if ((R >= ((4 * _i1) + 60)))
    {
      #pragma ivdep
      for (int  _i2 = 3; (_i2 < ((C / 8) - 5)); _i2 = (_i2 + 1))
      {
        Ux_outGPyramid_L2[((((_i1 + 1) - 7) * (((((C / 8) - 8) + 2) - 3) + 1)) + (_i2 - 3))] = ((0.25f * outGPyramid_L3[((((((_i1 + 1) / 2) - 1) - 3) * (((((C / 8) - 8) + 2) - 3) + 1)) + (-3 + _i2))]) + (0.75f * outGPyramid_L3[(((((_i1 + 1) / 2) - 3) * (((((C / 8) - 8) + 2) - 3) + 1)) + (-3 + _i2))]));
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 7; (_i1 < ((R / 4) - 13)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 7; (_i2 < ((C / 4) - 13)); _i2 = (_i2 + 1))
    {
      int  _ct90 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct91 = 0;
      int  _ct92 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct90: _ct91);
      int  _ct93 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct94 = 0;
      int  _ct95 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct93: _ct94);
      int  _ct96 = _ct95;
      int  _ct97 = 2;
      int  _ct98 = ((_ct92 < 2)? _ct96: _ct97);
      int  _ct99 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct100 = 0;
      int  _ct101 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct99: _ct100);
      int  _ct102 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct103 = 0;
      int  _ct104 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct102: _ct103);
      int  _ct105 = _ct104;
      int  _ct106 = 2;
      int  _ct107 = ((_ct101 < 2)? _ct105: _ct106);
      int  _ct108 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct109 = 0;
      int  _ct110 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct108: _ct109);
      int  _ct111 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct112 = 0;
      int  _ct113 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct111: _ct112);
      int  _ct114 = _ct113;
      int  _ct115 = 2;
      int  _ct116 = ((_ct110 < 2)? _ct114: _ct115);
      int  _ct117 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct118 = 0;
      int  _ct119 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct117: _ct118);
      int  _ct120 = (int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct121 = 0;
      int  _ct122 = (((int ) ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct120: _ct121);
      int  _ct123 = _ct122;
      int  _ct124 = 2;
      int  _ct125 = ((_ct119 < 2)? _ct123: _ct124);
      outLPyramid_L2[(((_i1 - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (_i2 - 7))] = (((1.0f - ((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct98)) * U_lPyramid_L2[(((_ct107 * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((-7 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1))) + (-7 + _i2))]) + (((D_inGPyramid_L2[(((-1 + _i1) * ((((C / 4) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct116) * U_lPyramid_L2[((((_ct125 + 1) * ((((((R / 4) - 16) + 2) - 7) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((-7 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1))) + (-7 + _i2))]));
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 15; (_i1 < ((R / 2) - 29)); _i1 = (_i1 + 2))
    {
      #pragma ivdep
      for (int  _i2 = 7; (_i2 < ((C / 4) - 13)); _i2 = (_i2 + 1))
      {
        Ux_lPyramid_L1[(((_i0 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((_i1 - 15) * (((((C / 4) - 16) + 2) - 7) + 1))) + (_i2 - 7))] = ((0.25f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))]));
      }
      if ((R >= ((2 * _i1) + 62)))
      {
        #pragma ivdep
        for (int  _i2 = 7; (_i2 < ((C / 4) - 13)); _i2 = (_i2 + 1))
        {
          Ux_lPyramid_L1[(((_i0 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + (((_i1 + 1) - 15) * (((((C / 4) - 16) + 2) - 7) + 1))) + (_i2 - 7))] = ((0.25f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + (((((_i1 + 1) / 2) - 1) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L2[(((_i0 * (((((R / 4) - 2) - 1) + 1) * ((((C / 4) - 2) - 1) + 1))) + ((((_i1 + 1) / 2) - 1) * ((((C / 4) - 2) - 1) + 1))) + (-1 + _i2))]));
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 7; (_i1 < ((R / 4) - 13)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 7; (_i2 < ((C / 4) - 13)); _i2 = (_i2 + 1))
    {
      float  _ct126 = ((0.25f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + (((_i2 / 2) - 1) - 3))]) + (0.75f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + ((_i2 / 2) - 3))]));
      float  _ct127 = ((0.25f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + (((_i2 / 2) + 1) - 3))]) + (0.75f * Ux_outGPyramid_L2[(((-7 + _i1) * (((((C / 8) - 8) + 2) - 3) + 1)) + ((_i2 / 2) - 3))]));
      float  _ct128 = (((_i2 % 2) == 0)? _ct126: _ct127);
      outGPyramid_L2[(((_i1 - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (_i2 - 7))] = (outLPyramid_L2[(((-7 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + _i2))] + _ct128);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 15; (_i1 < ((R / 2) - 29)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 15; (_i2 < ((C / 2) - 29)); _i2 = (_i2 + 1))
      {
        float  _ct129 = ((0.25f * Ux_lPyramid_L1[(((_i0 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1))) + (((_i2 / 2) - 1) - 7))]) + (0.75f * Ux_lPyramid_L1[(((_i0 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((_i2 / 2) - 7))]));
        float  _ct130 = ((0.25f * Ux_lPyramid_L1[(((_i0 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1))) + (((_i2 / 2) + 1) - 7))]) + (0.75f * Ux_lPyramid_L1[(((_i0 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1))) + ((_i2 / 2) - 7))]));
        float  _ct131 = (((_i2 % 2) == 0)? _ct129: _ct130);
        U_lPyramid_L1[(((_i0 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((_i1 - 15) * (((((C / 2) - 32) + 2) - 15) + 1))) + (_i2 - 15))] = (D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((-1 + _i1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))] - _ct131);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 15; (_i1 < ((R / 2) - 29)); _i1 = (_i1 + 2))
  {
    #pragma ivdep
    for (int  _i2 = 7; (_i2 < ((C / 4) - 13)); _i2 = (_i2 + 1))
    {
      Ux_outGPyramid_L1[(((_i1 - 15) * (((((C / 4) - 16) + 2) - 7) + 1)) + (_i2 - 7))] = ((0.25f * outGPyramid_L2[(((((_i1 / 2) + 1) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + _i2))]) + (0.75f * outGPyramid_L2[((((_i1 / 2) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + _i2))]));
    }
    if ((R >= ((2 * _i1) + 62)))
    {
      #pragma ivdep
      for (int  _i2 = 7; (_i2 < ((C / 4) - 13)); _i2 = (_i2 + 1))
      {
        Ux_outGPyramid_L1[((((_i1 + 1) - 15) * (((((C / 4) - 16) + 2) - 7) + 1)) + (_i2 - 7))] = ((0.25f * outGPyramid_L2[((((((_i1 + 1) / 2) - 1) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + _i2))]) + (0.75f * outGPyramid_L2[(((((_i1 + 1) / 2) - 7) * (((((C / 4) - 16) + 2) - 7) + 1)) + (-7 + _i2))]));
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 15; (_i1 < ((R / 2) - 29)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 15; (_i2 < ((C / 2) - 29)); _i2 = (_i2 + 1))
    {
      int  _ct132 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct133 = 0;
      int  _ct134 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct132: _ct133);
      int  _ct135 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct136 = 0;
      int  _ct137 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct135: _ct136);
      int  _ct138 = _ct137;
      int  _ct139 = 2;
      int  _ct140 = ((_ct134 < 2)? _ct138: _ct139);
      int  _ct141 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct142 = 0;
      int  _ct143 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct141: _ct142);
      int  _ct144 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct145 = 0;
      int  _ct146 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct144: _ct145);
      int  _ct147 = _ct146;
      int  _ct148 = 2;
      int  _ct149 = ((_ct143 < 2)? _ct147: _ct148);
      int  _ct150 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct151 = 0;
      int  _ct152 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct150: _ct151);
      int  _ct153 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct154 = 0;
      int  _ct155 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct153: _ct154);
      int  _ct156 = _ct155;
      int  _ct157 = 2;
      int  _ct158 = ((_ct152 < 2)? _ct156: _ct157);
      int  _ct159 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct160 = 0;
      int  _ct161 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct159: _ct160);
      int  _ct162 = (int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)));
      int  _ct163 = 0;
      int  _ct164 = (((int ) ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3))) > 0)? _ct162: _ct163);
      int  _ct165 = _ct164;
      int  _ct166 = 2;
      int  _ct167 = ((_ct161 < 2)? _ct165: _ct166);
      outLPyramid_L1[(((_i1 - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (_i2 - 15))] = (((1.0f - ((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct140)) * U_lPyramid_L1[(((_ct149 * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((-15 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1))) + (-15 + _i2))]) + (((D_inGPyramid_L1[(((-1 + _i1) * ((((C / 2) - 2) - 1) + 1)) + (-1 + _i2))] * (float ) (3)) - _ct158) * U_lPyramid_L1[((((_ct167 + 1) * ((((((R / 2) - 32) + 2) - 15) + 1) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((-15 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1))) + (-15 + _i2))]));
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 31; (_i1 < (R - 61)); _i1 = (_i1 + 2))
    {
      #pragma ivdep
      for (int  _i2 = 15; (_i2 < ((C / 2) - 29)); _i2 = (_i2 + 1))
      {
        Ux_lPyramid_L0[(((_i0 * ((-92 + R) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((_i1 - 31) * (((((C / 2) - 32) + 2) - 15) + 1))) + (_i2 - 15))] = ((0.25f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((((_i1 / 2) + 1) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + (((_i1 / 2) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))]));
      }
      if ((R >= (_i1 + 63)))
      {
        #pragma ivdep
        for (int  _i2 = 15; (_i2 < ((C / 2) - 29)); _i2 = (_i2 + 1))
        {
          Ux_lPyramid_L0[(((_i0 * ((-92 + R) * (((((C / 2) - 32) + 2) - 15) + 1))) + (((_i1 + 1) - 31) * (((((C / 2) - 32) + 2) - 15) + 1))) + (_i2 - 15))] = ((0.25f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + (((((_i1 + 1) / 2) - 1) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))]) + (0.75f * D_gPyramid_L1[(((_i0 * (((((R / 2) - 2) - 1) + 1) * ((((C / 2) - 2) - 1) + 1))) + ((((_i1 + 1) / 2) - 1) * ((((C / 2) - 2) - 1) + 1))) + (-1 + _i2))]));
        }
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 15; (_i1 < ((R / 2) - 29)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 15; (_i2 < ((C / 2) - 29)); _i2 = (_i2 + 1))
    {
      float  _ct168 = ((0.25f * Ux_outGPyramid_L1[(((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1)) + (((_i2 / 2) - 1) - 7))]) + (0.75f * Ux_outGPyramid_L1[(((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1)) + ((_i2 / 2) - 7))]));
      float  _ct169 = ((0.25f * Ux_outGPyramid_L1[(((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1)) + (((_i2 / 2) + 1) - 7))]) + (0.75f * Ux_outGPyramid_L1[(((-15 + _i1) * (((((C / 4) - 16) + 2) - 7) + 1)) + ((_i2 / 2) - 7))]));
      float  _ct170 = (((_i2 % 2) == 0)? _ct168: _ct169);
      outGPyramid_L1[(((_i1 - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (_i2 - 15))] = (outLPyramid_L1[(((-15 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + _i2))] + _ct170);
    }
  }
  for (int  _i0 = 0; (_i0 <= 3); _i0 = (_i0 + 1))
  {
    #pragma omp parallel for schedule(static)
    for (int  _i1 = 31; (_i1 < (R - 61)); _i1 = (_i1 + 1))
    {
      #pragma ivdep
      for (int  _i2 = 31; (_i2 < (C - 61)); _i2 = (_i2 + 1))
      {
        float  _ct171 = ((0.25f * Ux_lPyramid_L0[(((_i0 * ((-92 + R) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1))) + (((_i2 / 2) - 1) - 15))]) + (0.75f * Ux_lPyramid_L0[(((_i0 * ((-92 + R) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((_i2 / 2) - 15))]));
        float  _ct172 = ((0.25f * Ux_lPyramid_L0[(((_i0 * ((-92 + R) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1))) + (((_i2 / 2) + 1) - 15))]) + (0.75f * Ux_lPyramid_L0[(((_i0 * ((-92 + R) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1))) + ((_i2 / 2) - 15))]));
        float  _ct173 = (((_i2 % 2) == 0)? _ct171: _ct172);
        U_lPyramid_L0[(((_i0 * ((-92 + R) * (-92 + C))) + ((_i1 - 31) * (-92 + C))) + (_i2 - 31))] = (gPyramid_L0[(((_i0 * (R * C)) + (_i1 * C)) + _i2)] - _ct173);
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 31; (_i1 < (R - 61)); _i1 = (_i1 + 2))
  {
    #pragma ivdep
    for (int  _i2 = 15; (_i2 < ((C / 2) - 29)); _i2 = (_i2 + 1))
    {
      Ux_result_ref_gray[(((_i1 - 31) * (((((C / 2) - 32) + 2) - 15) + 1)) + (_i2 - 15))] = ((0.25f * outGPyramid_L1[(((((_i1 / 2) + 1) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + _i2))]) + (0.75f * outGPyramid_L1[((((_i1 / 2) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + _i2))]));
    }
    if ((R >= (_i1 + 63)))
    {
      #pragma ivdep
      for (int  _i2 = 15; (_i2 < ((C / 2) - 29)); _i2 = (_i2 + 1))
      {
        Ux_result_ref_gray[((((_i1 + 1) - 31) * (((((C / 2) - 32) + 2) - 15) + 1)) + (_i2 - 15))] = ((0.25f * outGPyramid_L1[((((((_i1 + 1) / 2) - 1) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + _i2))]) + (0.75f * outGPyramid_L1[(((((_i1 + 1) / 2) - 15) * (((((C / 2) - 32) + 2) - 15) + 1)) + (-15 + _i2))]));
      }
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 31; (_i1 < (R - 61)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 31; (_i2 < (C - 61)); _i2 = (_i2 + 1))
    {
      int  _ct174 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct175 = 0;
      int  _ct176 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct174: _ct175);
      int  _ct177 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct178 = 0;
      int  _ct179 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct177: _ct178);
      int  _ct180 = _ct179;
      int  _ct181 = 2;
      int  _ct182 = ((_ct176 < 2)? _ct180: _ct181);
      int  _ct183 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct184 = 0;
      int  _ct185 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct183: _ct184);
      int  _ct186 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct187 = 0;
      int  _ct188 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct186: _ct187);
      int  _ct189 = _ct188;
      int  _ct190 = 2;
      int  _ct191 = ((_ct185 < 2)? _ct189: _ct190);
      int  _ct192 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct193 = 0;
      int  _ct194 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct192: _ct193);
      int  _ct195 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct196 = 0;
      int  _ct197 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct195: _ct196);
      int  _ct198 = _ct197;
      int  _ct199 = 2;
      int  _ct200 = ((_ct194 < 2)? _ct198: _ct199);
      int  _ct201 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct202 = 0;
      int  _ct203 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct201: _ct202);
      int  _ct204 = (int ) ((img[((_i1 * C) + _i2)] * (float ) (3)));
      int  _ct205 = 0;
      int  _ct206 = (((int ) ((img[((_i1 * C) + _i2)] * (float ) (3))) > 0)? _ct204: _ct205);
      int  _ct207 = _ct206;
      int  _ct208 = 2;
      int  _ct209 = ((_ct203 < 2)? _ct207: _ct208);
      outLPyramid_L0[(((_i1 - 31) * (-92 + C)) + (_i2 - 31))] = (((1.0f - ((img[((_i1 * C) + _i2)] * (float ) (3)) - _ct182)) * U_lPyramid_L0[(((_ct191 * ((-92 + R) * (-92 + C))) + ((-31 + _i1) * (-92 + C))) + (-31 + _i2))]) + (((img[((_i1 * C) + _i2)] * (float ) (3)) - _ct200) * U_lPyramid_L0[((((_ct209 + 1) * ((-92 + R) * (-92 + C))) + ((-31 + _i1) * (-92 + C))) + (-31 + _i2))]));
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 31; (_i1 < (R - 61)); _i1 = (_i1 + 1))
  {
    #pragma ivdep
    for (int  _i2 = 31; (_i2 < (C - 61)); _i2 = (_i2 + 1))
    {
      float  _ct210 = ((0.25f * Ux_result_ref_gray[(((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1)) + (((_i2 / 2) - 1) - 15))]) + (0.75f * Ux_result_ref_gray[(((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1)) + ((_i2 / 2) - 15))]));
      float  _ct211 = ((0.25f * Ux_result_ref_gray[(((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1)) + (((_i2 / 2) + 1) - 15))]) + (0.75f * Ux_result_ref_gray[(((-31 + _i1) * (((((C / 2) - 32) + 2) - 15) + 1)) + ((_i2 / 2) - 15))]));
      float  _ct212 = (((_i2 % 2) == 0)? _ct210: _ct211);
      result_ref_gray[(((_i1 - 31) * (-92 + C)) + (_i2 - 31))] = (outLPyramid_L0[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] + _ct212);
    }
  }
  #pragma omp parallel for schedule(static)
  for (int  _i1 = 31; (_i1 < (R - 61)); _i1 = (_i1 + 1))
  {
    for (int  _i2 = 31; (_i2 < (C - 61)); _i2 = (_i2 + 1))
    {
      #pragma ivdep
      for (int  _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
      {
        float  _ct20426 = ((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f));
        float  _ct20427 = 0.0f;
        float  _ct20428 = ((((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f)) > 0.0f)? _ct20426: _ct20427);
        float  _ct20429 = ((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f));
        float  _ct20430 = 0.0f;
        float  _ct20431 = ((((result_ref_gray[(((-31 + _i1) * (-92 + C)) + (-31 + _i2))] * ((img_colour[(((_i1 * C * 3)) + _i2 * 3 + _i0)] / 65535.0f) + 0.01f)) / (img[((_i1 * C) + _i2)] + 0.01f)) > 0.0f)? _ct20429: _ct20430);
        float  _ct20432 = _ct20431;
        float  _ct20433 = 1.0f;
        float  _ct20434 = ((_ct20428 < 1.0f)? _ct20432: _ct20433);
        laplacian[((((_i1 - 31) * (-92 + C) * 3)) + (_i2 - 31) * 3 + _i0)] = (short unsigned int ) ((_ct20434 * 65535.0f));
      }
    }
  }
  free(img);
  free(img_colour);
  free(remapLUT);
  free(Dx_inGPyramid_L1);
  free(gPyramid_L0);
  free(D_inGPyramid_L1);
  free(Dx_gPyramid_L1);
  free(D_gPyramid_L1);
  free(Dx_inGPyramid_L2);
  free(D_inGPyramid_L2);
  free(Dx_gPyramid_L2);
  free(Ux_lPyramid_L0);
  free(D_gPyramid_L2);
  free(Dx_inGPyramid_L3);
  free(U_lPyramid_L0);
  free(D_inGPyramid_L3);
  free(Dx_gPyramid_L3);
  free(Ux_lPyramid_L1);
  free(outLPyramid_L0);
  free(D_gPyramid_L3);
  free(Dx_inGPyramid_L4);
  free(U_lPyramid_L1);
  free(D_inGPyramid_L4);
  free(Dx_gPyramid_L4);
  free(Ux_lPyramid_L2);
  free(outLPyramid_L1);
  free(D_gPyramid_L4);
  free(U_lPyramid_L2);
  free(Ux_lPyramid_L3);
  free(outLPyramid_L2);
  free(outLPyramid_L4);
  free(U_lPyramid_L3);
  free(Ux_outGPyramid_L3);
  free(outLPyramid_L3);
  free(outGPyramid_L3);
  free(Ux_outGPyramid_L2);
  free(outGPyramid_L2);
  free(Ux_outGPyramid_L1);
  free(outGPyramid_L1);
  free(Ux_result_ref_gray);
  free(result_ref_gray);
}
