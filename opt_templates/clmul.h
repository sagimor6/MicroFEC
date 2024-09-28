
#ifndef __CLMUL_H__
#define __CLMUL_H__

#include "../fec_common.h"

#if defined(FEC_HAS_CLMUL32)

#if (defined(__x86_64__) || defined(__i386__)) && defined(__PCLMUL__) && defined(__SSE2__)
#define _poly_t u32x4
#define _POLY_CLMUL(poly1, poly2) ((_poly_t)_mm_clmulepi64_si128((__m128i)(poly1), (__m128i)(poly2), 0))
#define _POLY_2VAL(val1, val2) ((_poly_t){(val1), (val2), 0, 0})
#define _POLY_EXTRACT(poly_val, typ, idx) (((typ __attribute__((vector_size (16))))(poly_val))[idx])
#define _POLY_VEC_SHIFT(poly_val, shift) ((poly_val) >> shift)
#elif (defined(__aarch64__) || defined(__arm__)) && defined(__ARM_FEATURE_AES)
// TODO: FEAT_PMULL is needed in processor
#define _poly_t poly64_t
#define _POLY_CLMUL(poly1, poly2) ((_poly_t)vmull_p64((poly1), (poly2)))
#elif defined(__sparc__) && defined(__VIS) && __VIS >= 0x300
#define _poly_t uint64_t
#define _POLY_CLMUL(poly1, poly2) ((_poly_t)__builtin_vis_xmulx((int64_t)(poly1), (int64_t)(poly2)))
#elif defined(__riscv__) && (defined(__riscv_zbc) || defined(__riscv_zbkc))
#if __riscv_xlen == 64
#define _poly_t uint64_t
#define _POLY_CLMUL(poly1, poly2) ((_poly_t)__builtin_riscv_clmul_64((poly1), (poly2)))
#elif __riscv_xlen == 32
#define _poly_t uint32_t
#define _POLY_CLMUL(poly1, poly2) ((_poly_t)__builtin_riscv_clmul_32((poly1), (poly2)))
#endif
#endif


#ifndef _POLY_2VAL
#ifdef FEC_HAS_CLMUL64
#define _POLY_2VAL(val1, val2) ((_poly_t)((val1) |  (((_poly_t)(val2)) << 32)))
#else
#define _POLY_2VAL(val1, val2) ((_poly_t)(val1))
#endif
#endif

#define _POLY_1VAL(val) _POLY_2VAL(val, 0)

#ifndef _POLY_EXTRACT
#define _POLY_EXTRACT(poly_val, typ, idx) ((typ)((poly_val) >> (sizeof(typ)*8*idx)))
#endif

#ifndef _POLY_VEC_SHIFT
#define _POLY_VEC_SHIFT(poly_val, shift) (((poly_val) >> (shift)) & (~((1ULL << 32) - (1ULL (32 - (shift))))))
#endif

#endif

#endif // __CLMUL_H__
