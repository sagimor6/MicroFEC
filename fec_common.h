#ifndef __FEC_COMMON_H__
#define __FEC_COMMON_H__

#include <stdint.h>

#include "micro_fec.h"

#if defined(PERF_DEBUG)
#define PERF_DEBUG_ATTRS __attribute__((noinline)) __attribute__((aligned(16)))
#else
#define PERF_DEBUG_ATTRS 
#endif

#ifndef MAX
#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif

#ifndef MIN
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

#define CALC_EXTRA_FOR_ALIGNED_MALLOC(typ) (__alignof__(typ) > __alignof__(max_align_t) ? __alignof__(typ) - __alignof__(max_align_t) : 0)

#define POLY_G ((fec_int_t)0b0000000000101011)

#if defined(FEC_HAS_64_INT_VEC) || defined (FEC_HAS_CLMUL32)

typedef uint8_t u8x16 __attribute__ ((vector_size (16)));
typedef uint16_t u16x8 __attribute__ ((vector_size (16)));
typedef uint32_t  u32x4 __attribute__ ((vector_size (16)));
typedef uint64_t  u64x2 __attribute__ ((vector_size (16)));

typedef uint8_t u8x32 __attribute__ ((vector_size (32)));
typedef uint16_t u16x16 __attribute__ ((vector_size (32)));
typedef uint32_t  u32x8 __attribute__ ((vector_size (32)));
typedef uint64_t  u64x4 __attribute__ ((vector_size (32)));

typedef int16_t i16x16 __attribute__ ((vector_size (32)));
typedef int16_t i16x8 __attribute__ ((vector_size (16)));

typedef uint8_t u8x8 __attribute__ ((vector_size (8)));
typedef uint16_t u16x4 __attribute__ ((vector_size (8)));
typedef uint32_t  u32x2 __attribute__ ((vector_size (8)));
typedef uint64_t  u64x1 __attribute__ ((vector_size (8)));

#if (defined(__x86_64__) || defined (__i386__))
#include <immintrin.h>
#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#endif

#if defined(__AVX2__)
#define my_mm256_extracti128_si256(x, y) _mm256_extracti128_si256((__m256i)x, y)
#elif defined(__AVX__)
#define my_mm256_extracti128_si256(x, y) _mm256_extractf128_si256((__m256i)x, y)
#else
#define my_mm256_extracti128_si256(x, y) ((u64x2){((u64x4)x)[y*2+0], ((u64x4)x)[y*2+1]})
#endif

#if defined(__SSE2__)
#define my_mm_shuffle_epi32(x, a, b, c, d) ((u32x4)_mm_shuffle_epi32((__m128i)x, _MM_SHUFFLE(d, c, b, a)))
#elif defined(__SSE__)
#define my_mm_shuffle_epi32(x, a, b, c, d) ((u32x4)_mm_shuffle_ps(x, x, _MM_SHUFFLE(d, c, b, a)))
#else
#define my_mm_shuffle_epi32(x, a, b, c, d) ((u32x4){((u32x4)x)[a], ((u32x4)x)[b], ((u32x4)x)[c], ((u32x4)x)[d]})
#endif

#endif


#endif // __FEC_COMMON_H__
