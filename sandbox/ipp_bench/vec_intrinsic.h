#include <immintrin.h>

using namespace std;

#define SSE 128
#define AVX 256
#define AVX_512 512

#ifndef VTYPE
#define VTYPE SSE
#endif

#if VTYPE == SSE
    #define VECTOR __m128
    #define VLEN 4
    #define ALIGN 16

    /* load, store */
    #define VECSPLAT(dst, src)      dst = _mm_load1_ps(src);
    #define VECSET(e3, e2, e1, e0)  _mm_set_ps(e3, e2, e1, e0);

    #define VECLOADU(src)			_mm_loadu_ps(src);
    #define VECLOAD(src)            _mm_load_ps(src);

    #define VECSTORE(src, dst)      _mm_store_ps(dst, src);
    #define VECSTOREU(src, dst)     _mm_storeu_ps(dst, src);

	#define VECSHUFFLE(a, b, imm)	_mm_shuffle_ps(a, b, imm);

    /* add, sub, mul, div */
    #define VECADD(vec1, vec2)      _mm_add_ps(vec1, vec2);
    #define VECSUB(vec1, vec2)      _mm_sub_ps(vec1, vec2);
    #define VECMUL(vec1, vec2)      _mm_mul_ps(vec1, vec2);

#elif VTYPE == AVX
    #define VECTOR __m256
    #define VLEN 8
    #define ALIGN 32

    /* load, store */
    #define VECSPLAT(dst, src)      dst = _mm256_loadu2_m128(src+4, src);
    #define VECSET(e7, e6, e5, e4, e3, e2, e1, e0) \
                                    _mm256_set_ps(e7, e6, e5, e4, e3, e2, e1, e0);

    #define VECLOADU(src)			_mm256_loadu_ps(src);
    #define VECLOAD(src)            _mm256_load_ps(src);

    #define VECSTORE(src, dst)      _mm256_store_ps(dst, src);
    #define VECSTOREU(src, dst)     _mm256_storeu_ps(dst, src);

	#define VECSHUFFLE(a, b, imm)	_mm256_shuffle_ps(a, b, imm);
    
    /* add, sub, mul, div */
    #define VECADD(vec1, vec2)      _mm256_add_ps(vec1, vec2);
    #define VECSUB(vec1, vec2)      _mm256_sub_ps(vec1, vec2);
    #define VECMUL(vec1, vec2)      _mm256_mul_ps(vec1, vec2);

#elif VTYPE == AVX_512
    #define VECTOR __m512
    #define VLEN 16
    #define ALIGN 64

    /* load, store */
    #define VECSPLAT(dst, src)      { \
                                    float __buf__[16]; \
                                    for(int i = 0; i < 16; i++) \
                                        __buf__[i] = *src; \
                                    dst = _mm512_load_ps(__buf__); \
                                    }
    #define VECSET(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0) \
                                    _mm512_set_ps(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0);

    #define VECLOADU(src)			_mm512_loadu_ps(src);
    #define VECLOAD(src)            _mm512_load_ps(src);

    #define VECSTORE(src, dst)      _mm512_store_ps(dst, src);
    #define VECSTOREU(src, dst)     _mm512_storeu_ps(dst, src);

	#define VECSHUFFLE(a, b, imm)	_mm512_shuffle_ps(a, b, imm);

    /* add, sub, mul, div */
    #define VECADD(vec1, vec2)      _mm512_add_ps(vec1, vec2);
    #define VECSUB(vec1, vec2)      _mm512_sub_ps(vec1, vec2);
    #define VECMUL(vec1, vec2)      _mm512_mul_ps(vec1, vec2);

#endif

