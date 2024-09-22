
#include <mmintrin.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#define PERF_DEBUG 
#ifdef PERF_DEBUG
#include <time.h>
#include <stdio.h>
#endif

#include "micro_fec.h"

#ifdef PERF_DEBUG
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

#define CALC_EXTRA_FOR_ALIGNED_MALLOC(typ) (__alignof__(typ) >= __alignof__(max_align_t) ? __alignof__(typ) - __alignof__(max_align_t) : 0)

#define POLY_G ((fec_int_t)0b0000000000101011)


// num must be > 0
static unsigned int log2_ceil(unsigned int num) {
    if (num == 1) {
        return 0;
    }
    return (sizeof(num)*8) - __builtin_clz(num - 1);
}

static fec_int_t poly_add(fec_int_t a, fec_int_t b) {
    return a ^ b;
}

#if defined(__x86_64__) || defined (__i386__)
#include <immintrin.h>

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

#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#endif

#if (defined(__PCLMUL__) && defined(__SSE2__)) || defined(_FEC_NO_OPT) || !(defined(__x86_64__) || defined(__i386__))
#define _FEC_USE_POLY_MUL
//#define _FEC_USE_POLY_MUL16
//#define _FEC_USE_POLY_MUL_CLMUL2
#elif defined(__AVX2__)
#define _FEC_USE_POLY_MUL16
#elif defined(__SSE2__)
#define _FEC_USE_POLY_MUL8
#elif defined(__x86_64__)
#define _FEC_USE_POLY_MUL4
#elif defined(__MMX__)
#define _FEC_USE_POLY_MUL4_MMX
#elif defined(__i386__)
#define _FEC_USE_POLY_MUL2
#else
#error all cases should be defined
#endif

#define ALIGN_UP(val, align) (((val) + (align) - 1) & (-((__typeof__(val))(align))))

#if defined(_FEC_USE_POLY_MUL) || defined(_FEC_USE_POLY_MUL_CLMUL2)
#define _FEC_ALIGN_SIZE_VAL sizeof(fec_int_t)
#elif defined(_FEC_USE_POLY_MUL16)
#define _FEC_ALIGN_SIZE_VAL 32
#elif defined(_FEC_USE_POLY_MUL8)
#define _FEC_ALIGN_SIZE_VAL 16
#elif defined(_FEC_USE_POLY_MUL4) || defined(_FEC_USE_POLY_MUL4_MMX)
#define _FEC_ALIGN_SIZE_VAL 8
#elif defined(_FEC_USE_POLY_MUL2)
#define _FEC_ALIGN_SIZE_VAL 4
#endif

#if (defined(__x86_64__) || defined(__i386__)) && defined(__PCLMUL__) && defined(__SSE2__)
#define _poly_t u32x4
#define _POLY_INT_TO_POLY(val) ((_poly_t){(val), 0, 0, 0})
#define _POLY_TO_FEC_INT(poly_val) (((u16x8)(poly_val))[0])
#define _POLY_SINGLE_MUL(poly1, poly2) ((_poly_t)_mm_clmulepi64_si128((__m128i)(poly1), (__m128i)(poly2), 0))
#define _POLY_TWO_INT_TO_POLY(val1, val2) ((_poly_t){(val1), (val2), 0, 0})
#define _POLY_GET_LOWER_32(poly_val) ((poly_val)[0])
#define _POLY_GET_UPPER_32(poly_val) ((poly_val)[1])
#define _POLY_GET_LOWER_FEC_INT(poly_val) (((u16x8)(poly_val))[0])
#define _POLY_GET_UPPER_FEC_INT(poly_val) (((u16x8)(poly_val))[2])
#define _POLY_TWO_SHIFT_RIGHT(poly_val, shift) ((poly_val) >> shift)
#elif (defined(__aarch64__) || defined(__arm__)) && defined(__ARM_FEATURE_AES)
// TODO: FEAT_PMULL is needed in processor
#define _poly_t poly64_t
#define _POLY_SINGLE_MUL(poly1, poly2) ((_poly_t)vmull_p64((poly1), (poly2)))
#elif defined(__sparc__) && defined(__VIS) && __VIS >= 0x300
#define _poly_t uint64_t
#define _POLY_SINGLE_MUL(poly1, poly2) ((_poly_t)__builtin_vis_xmulx((int64_t)(poly1), (int64_t)(poly2)))
#elif defined(__riscv__) && (defined(__riscv_zbc) || defined(__riscv_zbkc))
#if __riscv_xlen == 64
#define _poly_t uint64_t
#define _POLY_SINGLE_MUL(poly1, poly2) ((_poly_t)__builtin_riscv_clmul_64((poly1), (poly2)))
#elif __riscv_xlen == 32
#define _poly_t uint32_t
#define _POLY_SINGLE_MUL(poly1, poly2) ((_poly_t)__builtin_riscv_clmul_32((poly1), (poly2)))
#endif
#endif

#ifdef FEC_HAS_CLMUL32
#ifndef _POLY_INT_TO_POLY
#define _POLY_INT_TO_POLY(val) ((_poly_t)(val))
#endif
#ifndef _POLY_TO_FEC_INT
#define _POLY_TO_FEC_INT(poly_val) ((fec_int_t)(poly_val))
#endif
#endif

#ifdef FEC_HAS_CLMUL64
#ifndef _POLY_TWO_INT_TO_POLY
#define _POLY_TWO_INT_TO_POLY(val1, val2) ((_poly_t)((val1) |  (((_poly_t)(val2)) << 32)))
#endif
#ifndef _POLY_GET_LOWER_32
#define _POLY_GET_LOWER_32(poly_val) ((uint32_t)(poly_val))
#endif
#ifndef _POLY_GET_UPPER_32
#define _POLY_GET_UPPER_32(poly_val) ((uint32_t)((poly_val) >> 32))
#endif
#ifndef _POLY_GET_LOWER_FEC_INT
#define _POLY_GET_LOWER_FEC_INT(poly_val) ((fec_int_t)(poly_val))
#endif
#ifndef _POLY_GET_UPPER_FEC_INT
#define _POLY_GET_UPPER_FEC_INT(poly_val) ((fec_int_t)((poly_val) >> 32))
#endif
#ifndef _POLY_TWO_SHIFT_RIGHT
#define _POLY_TWO_SHIFT_RIGHT(poly_val, shift) (((poly_val) >> (shift)) & (~((1ULL << 32) - (1ULL (32 - (shift))))))
#endif
#endif


static fec_int_t PERF_DEBUG_ATTRS poly_mul(fec_int_t a, fec_int_t b) {
#if defined(FEC_HAS_CLMUL32)

    _poly_t _a = _POLY_INT_TO_POLY(a);
    _poly_t _b = _POLY_INT_TO_POLY(b);
    _poly_t _poly = _POLY_INT_TO_POLY(POLY_G);

    _poly_t _c;
    _c = _POLY_SINGLE_MUL(_a, _b);
    _poly_t _d;
    _d = _c >> 16;
    _d = _POLY_SINGLE_MUL(_d, _poly);
    _c ^= _d;
    _d >>= 16;
    _d = _POLY_SINGLE_MUL(_d, _poly);
    _c ^= _d;
    
    return _POLY_TO_FEC_INT(_c);
#else
    size_t i;
    fec_int_t res = 0;

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        res = poly_add((res << 1), POLY_G & (-(res >> (sizeof(fec_int_t)*8 - 1)))); // poly left shift 1
        res = poly_add(res, a & (-(b >> (sizeof(fec_int_t)*8 - 1)))); // add a if current bit in b is 1
        b <<= 1;
    }

    return res;
#endif
}

#if defined(_FEC_USE_POLY_MUL4)
static uint64_t PERF_DEBUG_ATTRS poly_mul4(uint64_t a, uint64_t b) {
    size_t i;
    uint64_t res = 0;
    uint64_t shift_mask = (1ULL<<16) | (1ULL<<32) | (1ULL<<48);
    const uint64_t msb_mask = (1ULL<<(15+16*0)) | (1ULL<<(15+16*1)) | (1ULL<<(15+16*2)) | (1ULL<<(15+16*3));
    const uint64_t poly_g = (msb_mask >> 15) * POLY_G;

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        uint64_t shifted_res = (res << 1) & (~shift_mask);
        uint64_t res_msbs = (res & msb_mask);
        res = shifted_res ^ (((res_msbs << 1) - (res_msbs >> 15)) & poly_g); // poly left shift 1
        uint64_t b_msbs = b & msb_mask;
        res ^= (((b_msbs << 1) - (b_msbs >> 15)) & a); // add a if current bit in b is 1
        b <<= 1;
        b &= (~shift_mask);
    }

    return res;
}
#endif

#if defined(_FEC_USE_POLY_MUL2)
static uint32_t PERF_DEBUG_ATTRS poly_mul2_2(uint32_t a, uint32_t b) {
    size_t i;
    uint32_t res = 0;
    uint32_t shift_mask = (1<<16);
    const uint32_t msb_mask = (1ULL<<(15+16*0)) | (1ULL<<(15+16*1));
    const uint32_t poly_g = (msb_mask >> 15) * POLY_G;

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        uint32_t shifted_res = (res << 1) & (~shift_mask);
        uint32_t res_msbs = (res & msb_mask);
        res = shifted_res ^ (((res_msbs << 1) - (res_msbs >> 15)) & poly_g); // poly left shift 1
        uint32_t b_msbs = b & msb_mask;
        res ^= (((b_msbs << 1) - (b_msbs >> 15)) & a); // add a if current bit in b is 1
        b <<= 1;
        b &= (~shift_mask);
    }

    return res;
}
#endif

#if defined(_FEC_USE_POLY_MUL4_MMX)
static __m64 PERF_DEBUG_ATTRS poly_mul4_mmx(__m64 a, __m64 b) {
    size_t i;
    __m64 res = {0};
    const __m64 poly_g = (__m64)((u16x4){POLY_G,POLY_G,POLY_G,POLY_G});

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        res = _mm_xor_si64(_mm_slli_pi16(res, 1), _mm_and_si64(poly_g, _mm_srai_pi16(res, 15))); // poly left shift 1
        res = _mm_xor_si64(res, _mm_and_si64(a, _mm_srai_pi16(b, 15))); // add a if current bit in b is 1
        b = _mm_slli_pi16(b, 1);
    }

    return res;
}
#endif

#if defined(_FEC_USE_POLY_MUL8)
static u16x8 PERF_DEBUG_ATTRS poly_mul8(u16x8 a, u16x8 b) {
    size_t i;
    u16x8 res = {0};
    u16x8 poly_g = {POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G};

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        res = (res << 1) ^ (poly_g & (-(res >> 15))); // poly left shift 1
        res ^= a & (-(b >> 15)); // add a if current bit in b is 1
        b <<= 1;
    }

    return res;
}
#endif

#if defined(_FEC_USE_POLY_MUL16)
static u16x16 PERF_DEBUG_ATTRS poly_mul16(u16x16 a, u16x16 b) {
    size_t i;
    u16x16 res = {0};
    u16x16 poly_g = {POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,
                    POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G,POLY_G};

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        res = (res << 1) ^ (poly_g & (-(res >> 15))); // poly left shift 1
        res ^= a & (-(b >> 15)); // add a if current bit in b is 1
        b <<= 1;
    }

    return res;
}
#endif

// TODO: this doesn't seem to increase performance
#if defined(_FEC_USE_POLY_MUL_CLMUL2)
static uint32_t PERF_DEBUG_ATTRS poly_mul2(fec_int_t a1, fec_int_t a2, fec_int_t b1, fec_int_t b2) {
    u64x2 _a;
    _a[0] = a1 | ((uint64_t)a2<<32);
    u64x2 _b;
    _b[0] = b1 | ((uint64_t)b2<<32);
    u64x2 _poly;
    _poly[0] = POLY_G;

    u16x8 _c;
    _c = (u16x8)_mm_clmulepi64_si128((__m128i)_a, (__m128i)_b, 0);
    _c = (u16x8)_mm_shuffle_epi32((__m128i)_c, _MM_SHUFFLE(3, 2, 2, 0));
    u32x4 _d;
    _d = (((u32x4)_c) >> 16);
    _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    _c ^= (u16x8)_d;
    _d >>= 16;
    _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    _c ^= (u16x8)_d;
    
    return _c[0] | ((uint32_t)_c[2]) << 16;

    // u64x2 _a;
    // _a[0] = a1;
    // a2 = a2;
    // u64x2 _b;
    // _b[0] = b1 | ((uint64_t)b2<<32);
    // u64x2 _poly;
    // _poly[0] = POLY_G;

    // u16x8 _c;
    // _c = (u16x8)_mm_clmulepi64_si128((__m128i)_a, (__m128i)_b, 0);
    // u32x4 _d;
    // _d = (((u32x4)_c) >> 16);
    // _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    // _c ^= (u16x8)_d;
    // _d >>= 16;
    // _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    // _c ^= (u16x8)_d;
    
    // return _c[0] | ((uint32_t)_c[2]) << 16;
}
#endif

static fec_int_t poly_pow(fec_int_t a, fec_int_t n) {
    fec_int_t res = 1;
    fec_int_t a_2_i = a;
    while (n != 0) {
        if ((n & 1) != 0) {
            res = poly_mul(res, a_2_i);
        }
        a_2_i = poly_mul(a_2_i, a_2_i);
        n >>= 1;
    }

    return res;
}

static fec_int_t poly_inv(fec_int_t a) {
    fec_int_t res = poly_pow(a, (fec_int_t)-2);
    return res;
}

static inline fec_int_t _fec_inv(const fec_inv_cache_t *inv_cache, fec_int_t a) {
    return inv_cache->inv_arr[a - 1];
}

fec_status_t fec_inv_cache_init(fec_inv_cache_t *inv_cache, fec_idx_t n, fec_idx_t k) {
    fec_idx_t n_k_1 = 0;

    if(n != 0 || k != 0) {
        n_k_1 = n + k - 1;
        if (n_k_1 < n) {
            return FEC_STATUS_INVALID_PARAMS;
        }
    }

    return fec_inv_cache_init_raw(inv_cache, n_k_1);
}

fec_status_t fec_inv_cache_init_raw(fec_inv_cache_t *inv_cache, fec_idx_t n_k_1) {

    if (n_k_1 > (((fec_idx_t)1)<<(sizeof(fec_int_t)*8))) {
        return FEC_STATUS_INVALID_PARAMS;
    }

    inv_cache->inv_arr_sz = n_k_1;

    if (n_k_1 > 1) {
        inv_cache->inv_arr_sz = (1 << log2_ceil(n_k_1));
        size_t arr_size = inv_cache->inv_arr_sz - 1;
        inv_cache->inv_arr = malloc(sizeof(inv_cache->inv_arr[0]) * arr_size);
        if (inv_cache->inv_arr == NULL) {
            return FEC_STATUS_OUT_OF_MEMORY;
        }

        fec_int_t i;
        for (i = 0; i < arr_size; i++) {
            inv_cache->inv_arr[i] = poly_inv(i + 1);
        }

    } else {
        inv_cache->inv_arr_sz = 1;
        inv_cache->inv_arr = NULL;
    }

    return FEC_STATUS_SUCCESS;
}

void fec_inv_cache_destroy(fec_inv_cache_t *inv_cache) {
    if (inv_cache->inv_arr != NULL) {
        free(inv_cache->inv_arr);
    }
}

fec_status_t fec_tx_init(fec_tx_state_t *tx_state, fec_idx_t n, size_t pak_len) {
    memset(tx_state, 0, sizeof(*tx_state));

    if(n > (((fec_idx_t)1)<<(sizeof(fec_int_t)*8))) {
        return FEC_STATUS_INVALID_PARAMS;
    }

    tx_state->n = n;
    tx_state->pak_len = pak_len;
    tx_state->paks = calloc(n, sizeof(tx_state->paks[0]));
    if (tx_state->paks == NULL) {
        fec_tx_destroy(tx_state);
        return FEC_STATUS_OUT_OF_MEMORY;
    }

#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_TX_OPT)
    tx_state->tmp_pak = malloc(pak_len * sizeof(tx_state->tmp_pak[0]) + CALC_EXTRA_FOR_ALIGNED_MALLOC(tx_state->tmp_pak[0]));
    if (tx_state->tmp_pak == NULL) {
        fec_tx_destroy(tx_state);
        return FEC_STATUS_OUT_OF_MEMORY;
    }
#endif

    return FEC_STATUS_SUCCESS;
}

void fec_tx_destroy(fec_tx_state_t *tx_state) {
    if (tx_state->paks != NULL) {
        free(tx_state->paks);
    }
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_TX_OPT)
    if (tx_state->tmp_pak != NULL) {
        free(tx_state->tmp_pak);
    }
#endif
}


#ifndef FEC_USER_GIVEN_BUFFER
fec_status_t fec_rx_init(fec_rx_state_t *rx_state, fec_idx_t n, fec_idx_t k, size_t pak_len)
#else
fec_status_t fec_rx_init(fec_rx_state_t *rx_state, fec_idx_t n, fec_idx_t k, size_t pak_len, void* dest_buf)
#endif
{
    memset(rx_state, 0, sizeof(*rx_state));

    if (n == 0) {
        return FEC_STATUS_SUCCESS;
    }

    fec_idx_t n_k_1 = n + k - 1;
    if (n_k_1 < k || n_k_1 > (((fec_idx_t)1)<<(sizeof(fec_int_t)*8))) {
        return FEC_STATUS_INVALID_PARAMS;
    }
    // if k == 0 then all k-1 are in MIN so the underflow will do nothing

#define MALLOC_ATTR(name, size) \
    rx_state->name = malloc(sizeof(rx_state->name[0]) * (size) + CALC_EXTRA_FOR_ALIGNED_MALLOC(rx_state->name[0])); \
    if (rx_state->name == NULL) { \
        fec_rx_destroy(rx_state); \
        return FEC_STATUS_OUT_OF_MEMORY; \
    }
#define CALLOC_ATTR(name, size) \
    rx_state->name = calloc((size), sizeof(rx_state->name[0])); \
    if (rx_state->name == NULL) { \
        fec_rx_destroy(rx_state); \
        return FEC_STATUS_OUT_OF_MEMORY; \
    }

    rx_state->n = n;
    rx_state->k = k;
    rx_state->pak_len = pak_len;

#ifndef FEC_LARGE_K
    CALLOC_ATTR(info_paks, n);
    CALLOC_ATTR(redundancy_paks, k);
    MALLOC_ATTR(present_x, MIN(k - 1, n));
#else
#ifndef FEC_USER_GIVEN_BUFFER
    MALLOC_ATTR(pak_arr, n);
#else
    rx_state->pak_buffer = (unaligend_fec_int_t*)dest_buf;
#endif

    CALLOC_ATTR(received_paks_bitmap, (n + k + (8-1))/8);
    MALLOC_ATTR(pak_xy_arr, n);
#endif

    MALLOC_ATTR(missing_y, MIN(k, n));

#ifdef FEC_MIN_MEM
    MALLOC_ATTR(tmp_recovered_ints, MIN(k, n));
#else
    MALLOC_ATTR(pi_xy_div_xx, MIN(k - 1, n));
    MALLOC_ATTR(pi_yx_div_yy, MIN(k, n));
    MALLOC_ATTR(tmp_vec_redundancy, MIN(k - 1, n));
    MALLOC_ATTR(tmp_vec_info, n - 1);
    MALLOC_ATTR(present_y, n - 1);
    MALLOC_ATTR(pi_ycomp_y_div_ycomp_x, n - 1);
#endif

#undef MALLOC_ATTR
#undef CALLOC_ATTR

    rx_state->num_info = 0;
    rx_state->num_redundant = 0;
    rx_state->max_x = n - 1;
#if defined(FEC_LARGE_K)
#ifndef FEC_USER_GIVEN_BUFFER
    rx_state->ones_pak = NULL;
#else
    rx_state->has_one_pak = false;
#endif
#endif

    return FEC_STATUS_SUCCESS;
}

void fec_rx_reset(fec_rx_state_t *rx_state) {
    if (rx_state->n == 0) {
        return;
    }
#ifndef FEC_LARGE_K
    memset(rx_state->info_paks, 0, rx_state->n * sizeof(rx_state->info_paks[0]));
    memset(rx_state->redundancy_paks, 0, rx_state->k * sizeof(rx_state->redundancy_paks[0]));
#else
    memset(rx_state->received_paks_bitmap, 0, ((rx_state->n + rx_state->k + (8-1))/8) * sizeof(rx_state->received_paks_bitmap[0]));
#endif
    rx_state->num_info = 0;
    rx_state->num_redundant = 0;
    rx_state->max_x = rx_state->n - 1;
#ifdef FEC_LARGE_K
#ifndef FEC_USER_GIVEN_BUFFER
    rx_state->ones_pak = NULL;
#else
    rx_state->has_one_pak = false;
#endif
#endif
}

void fec_rx_destroy(fec_rx_state_t *rx_state) {
#define FREE_ATTR(name) \
    if (rx_state->name != NULL) { \
        free(rx_state->name); \
    }
    
#ifndef FEC_LARGE_K
    FREE_ATTR(info_paks);
    FREE_ATTR(redundancy_paks);
    FREE_ATTR(present_x);
#else
#ifndef FEC_USER_GIVEN_BUFFER
    FREE_ATTR(pak_arr);
#endif
    FREE_ATTR(received_paks_bitmap);
    FREE_ATTR(pak_xy_arr);
#endif

    FREE_ATTR(missing_y);

#ifdef FEC_MIN_MEM
    FREE_ATTR(tmp_recovered_ints);
#else
    FREE_ATTR(pi_xy_div_xx);
    FREE_ATTR(pi_yx_div_yy);
    FREE_ATTR(tmp_vec_info);
    FREE_ATTR(tmp_vec_redundancy);
    FREE_ATTR(present_y);
    FREE_ATTR(pi_ycomp_y_div_ycomp_x);
#endif
#undef FREE_ATTR
}

fec_status_t fec_tx_add_info_pak(fec_tx_state_t *tx_state, const void* pak, fec_idx_t idx) {
    if (idx >= tx_state->n) {
        return FEC_STATUS_INVALID_PARAMS;
    }
    tx_state->paks[idx] = (const unaligend_fec_int_t*)pak;
    return FEC_STATUS_SUCCESS;
}

static bool _fec_can_recover(const fec_rx_state_t *rx_state) {
#if defined(FEC_LARGE_K) && !defined(FEC_USER_GIVEN_BUFFER)
    return (rx_state->num_info + rx_state->num_redundant + (rx_state->ones_pak != NULL) >= rx_state->n);
#else
    return (rx_state->num_info + rx_state->num_redundant >= rx_state->n); // we count the ones pak in the redundant
#endif
}

#ifndef FEC_LARGE_K

fec_status_t fec_rx_is_pak_needed(fec_rx_state_t *rx_state, fec_idx_t idx) {
    fec_int_t n = rx_state->n;

    if (idx >= n + rx_state->k) {
        return FEC_STATUS_INVALID_PARAMS;
    }

    if (_fec_can_recover(rx_state)) {
        return FEC_STATUS_CAN_DROP_ALREADY_RECOVERABLE;
    }

    if (idx < n) {
        if (rx_state->info_paks[idx] != NULL) {
            return FEC_STATUS_CAN_DROP_DUP_PAK;
        }
    } else {
        if (rx_state->redundancy_paks[idx - n] != NULL) {
            return FEC_STATUS_CAN_DROP_DUP_PAK;
        }
    }

    return FEC_STATUS_SUCCESS;
}

fec_status_t fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_idx_t idx) {
    fec_int_t n = rx_state->n;
    fec_status_t status;

    status = fec_rx_is_pak_needed(rx_state, idx);
    if (status != FEC_STATUS_SUCCESS) {
        return status;
    }

    if (idx < n) {
        rx_state->info_paks[idx] = (unaligend_fec_int_t*)pak;
        rx_state->num_info++;
    } else {
        rx_state->redundancy_paks[idx - n] = (unaligend_fec_int_t*)pak;
        if (idx != n) {
            bool has_one_row = (rx_state->redundancy_paks[0] != NULL);
            rx_state->present_x[rx_state->num_redundant - has_one_row] = idx - 1;
        }
        rx_state->num_redundant++;
    }

    if (_fec_can_recover(rx_state)) {
        return FEC_STATUS_SUCCESS;
    } else {
        return FEC_STATUS_MORE_PACKETS_NEEDED;
    }
}

#else

fec_status_t fec_rx_is_pak_needed(fec_rx_state_t *rx_state, fec_idx_t idx) {
    fec_int_t n = rx_state->n;

    if (idx >= n + rx_state->k) {
        return FEC_STATUS_INVALID_PARAMS;
    }

    if (_fec_can_recover(rx_state)) {
        return FEC_STATUS_CAN_DROP_ALREADY_RECOVERABLE;
    }

    if ((rx_state->received_paks_bitmap[idx / 8] & (1<<(idx & (8-1)))) != 0) {
        return FEC_STATUS_CAN_DROP_DUP_PAK; // already received
    }

    return FEC_STATUS_SUCCESS;
}

fec_status_t fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_idx_t idx) {
    fec_int_t n = rx_state->n;
    fec_status_t status;

    status = fec_rx_is_pak_needed(rx_state, idx);
    if (status != FEC_STATUS_SUCCESS) {
        return status;
    }

    rx_state->received_paks_bitmap[idx / 8] |= (1<<(idx & (8-1)));

    if (idx < n) {
#ifndef FEC_USER_GIVEN_BUFFER
        rx_state->pak_arr[rx_state->num_info] = (unaligend_fec_int_t*)pak;
#else
        // we use memcpy and not memmove cause the user shouldn't be using the supplied packet buffer
        memcpy(&rx_state->pak_buffer[rx_state->num_info*rx_state->pak_len], pak, rx_state->pak_len * sizeof(fec_int_t));
#endif
        rx_state->pak_xy_arr[rx_state->num_info] = idx;
        rx_state->num_info++;
    } else if (idx == n) {
#ifndef FEC_USER_GIVEN_BUFFER
        rx_state->ones_pak = pak;
#else
        rx_state->ones_pak_idx = rx_state->num_redundant;
        memcpy(&rx_state->pak_buffer[(n - 1 - rx_state->num_redundant)*rx_state->pak_len], pak, rx_state->pak_len * sizeof(fec_int_t));
        rx_state->pak_xy_arr[n - 1 - rx_state->num_redundant] = 0; // invalid x
        rx_state->num_redundant++;
        rx_state->has_one_pak = true;
#endif
    } else {
#ifndef FEC_USER_GIVEN_BUFFER
        rx_state->pak_arr[n - 1 - rx_state->num_redundant] = (unaligend_fec_int_t*)pak;
#else
        memcpy(&rx_state->pak_buffer[(n - 1 - rx_state->num_redundant)*rx_state->pak_len], pak, rx_state->pak_len * sizeof(fec_int_t));
#endif
        rx_state->pak_xy_arr[n - 1 - rx_state->num_redundant] = idx - 1;
        rx_state->max_x = MAX(rx_state->max_x, (fec_int_t)(idx - 1));
        rx_state->num_redundant++;
    }

    if (_fec_can_recover(rx_state)) {
        return FEC_STATUS_SUCCESS;
    } else {
        return FEC_STATUS_MORE_PACKETS_NEEDED;
    }
}

#endif

#if defined(FEC_HAS_CLMUL64)
void fec_tx_init_perf_arr(uint32_t* restrict out_pak, size_t pak_len) {
    memset(out_pak, 0, pak_len*sizeof(out_pak[0]));
}

void PERF_DEBUG_ATTRS fec_tx_col_op(uint32_t* restrict out_pak, const unaligend_fec_int_t* restrict pak, size_t pak_len, fec_int_t a) {
    size_t j;

    _poly_t _a = _POLY_INT_TO_POLY(a);

    for (j = 0; j < pak_len - 1; j += 2) {
        _poly_t _b = _POLY_TWO_INT_TO_POLY(pak[j], pak[j+1]);

        _poly_t _c = _POLY_SINGLE_MUL(_a, _b);

        out_pak[j] ^= _POLY_GET_LOWER_32(_c);
        out_pak[j+1] ^= _POLY_GET_UPPER_32(_c);
    }
    
    if (j == pak_len - 1) {
        _poly_t _b = _POLY_INT_TO_POLY(pak[j]);

        _poly_t _c = _POLY_SINGLE_MUL(_a, _b);

        out_pak[j] ^= _POLY_GET_LOWER_32(_c);
    }
}

void fec_tx_col_perf_to_norm(unaligend_fec_int_t* restrict out_pak, const uint32_t* restrict perf_pak, size_t pak_len) {
    size_t j;
    _poly_t _poly = _POLY_INT_TO_POLY(POLY_G);

    for (j = 0; j < pak_len - 1; j += 2) {

        _poly_t _c = _POLY_TWO_INT_TO_POLY(perf_pak[j], perf_pak[j+1]);
        _poly_t _d = _POLY_TWO_SHIFT_RIGHT(_c, 16);
        _d = _POLY_SINGLE_MUL(_d, _poly);
        _c ^= _d;
        _d = _POLY_TWO_SHIFT_RIGHT(_d, 16);
        _d = _POLY_SINGLE_MUL(_d, _poly);
        _c ^= _d;
        
        out_pak[j] = _POLY_GET_LOWER_FEC_INT(_c);
        out_pak[j+1] = _POLY_GET_UPPER_FEC_INT(_c);
    }

    if (j == pak_len - 1) {
        _poly_t _c = _POLY_INT_TO_POLY(perf_pak[j]);
        _poly_t _d = _c >> 16;
        _d = _POLY_SINGLE_MUL(_d, _poly);
        _c ^= _d;
        _d >>= 16;
        _d = _POLY_SINGLE_MUL(_d, _poly);
        _c ^= _d;
        
        out_pak[j] = _POLY_TO_FEC_INT(_c);
    }
}
#endif

#if defined(__AVX2__)
#define my_mm256_extracti128_si256 _mm256_extracti128_si256
#elif defined(__AVX__)
#define my_mm256_extracti128_si256 _mm256_extractf128_si256
#else
#define my_mm256_extracti128_si256(x, y) ((u64x2){((u64x4)x)[y*2+0], ((u64x4)x)[y*2+1]})
#endif

#if defined(__SSE2__)
#define my_mm_shuffle_epi32(x, a, b, c, d) _mm_shuffle_epi32(x, _MM_SHUFFLE(d, c, b, a))
#elif defined(__SSE__)
#define my_mm_shuffle_epi32(x, a, b, c, d) ((u32x4)_mm_shuffle_ps(x, x, _MM_SHUFFLE(d, c, b, a)))
#else
#define my_mm_shuffle_epi32(x, a, b, c, d) ((u32x4){((u32x4)x)[a], ((u32x4)x)[b], ((u32x4)x)[c], ((u32x4)x)[d]})
#endif

void fec_tx_init_perf_arr_avx(u16x16* restrict out_pak, size_t pak_len) {
    memset(out_pak, 0, pak_len*sizeof(out_pak[0]));
}

#if defined(FEC_HAS_128_INT_VEC) || (defined(FEC_HAS_64_INT_VEC) && !(defined(__x86_64__) || defined(__i386__)))
void PERF_DEBUG_ATTRS fec_tx_col_op_avx(u16x16* restrict out_pak, const unaligend_fec_int_t* restrict pak, size_t pak_len, fec_int_t a) {
    fec_idx_t i;

    u16x16 shifts = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    u16x16 _a = (u16x16)(((i16x16)(a << shifts)) >> 15);

    for (i = 0; i < pak_len; i++) {
        u16x16 _c = _a & pak[i];
        out_pak[i] ^= _c;
    }
}
#endif

#if defined(FEC_HAS_128_INT_VEC)
void fec_tx_col_perf_to_norm_avx(unaligend_fec_int_t* restrict out_pak, const u16x16* restrict perf_pak, size_t pak_len) {
    fec_idx_t i;

    u32x8 shifts1 = {0, 2, 4, 6, 8, 10, 12, 14};
    u32x8 shifts2 = {1, 3, 5, 7, 9, 11, 13, 15};

    for (i = 0; i < pak_len; i++) {
        u16x16 res = perf_pak[i];
        u32x8 res1 = (((u32x8)res) & 0xFFFF) << shifts1;
        u32x8 res2 = (((u32x8)res) >> 16) << shifts2;

        u32x8 res3 = res1 ^ res2;
        u32x4 res4 = (u32x4)(my_mm256_extracti128_si256((__m256i)res3, 0) ^ my_mm256_extracti128_si256((__m256i)res3, 1));
        u32x4 res5 = res4 ^ my_mm_shuffle_epi32(res4, 1, 2, 3, 0);
        u32x4 res6 = res5 ^ my_mm_shuffle_epi32(res5, 2, 3, 0, 1); // res6 holds in each cell the int32 res

        u32x4 res7 = (res6 >> 16) << ((u32x4){0, 1, 3, 5}); // mul upper by POLY_G
        
        u32x4 res8 = res7 ^ my_mm_shuffle_epi32(res7, 1, 2, 3, 0);
        u32x4 res9 = res8 ^ my_mm_shuffle_epi32(res8, 2, 3, 0, 1);

        u32x4 res10 = (res9 >> 16) << ((u32x4){0, 1, 3, 5}); // mul upper by POLY_G
        
        u32x4 res11 = res10 ^ my_mm_shuffle_epi32(res10, 1, 2, 3, 0);
        u32x4 res12 = res11 ^ my_mm_shuffle_epi32(res11, 2, 3, 0, 1);

        out_pak[i] = ((u16x8)(res6 ^ res9 ^ res12))[0];

    }
}
#endif


#ifdef FEC_HAS_64_INT_VEC

void fec_tx_init_perf_arr_mmx(__m64* restrict out_pak, size_t pak_len) {
    memset(out_pak, 0, pak_len*sizeof(out_pak[0])*4);
}

#if (defined(__x86_64__) || defined(__i386__)) && defined(__MMX__)
void PERF_DEBUG_ATTRS __attribute__((visibility("default"))) __attribute__((used)) fec_tx_col_op_mmx(__m64* restrict out_pak, const unaligend_fec_int_t* restrict pak, size_t pak_len, fec_int_t a) {
    fec_idx_t i;

    int j;
    __m64 _a[4];
    const __m64 __a = _mm_set1_pi16(a);
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        _a[j] = _mm_srai_pi16(_mm_mullo_pi16(__a, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}), 15);
        //_a[j] = _mm_cmpgt_pi16((__m64){0}, _mm_mullo_pi16(__a, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}));
    }

    for (i = 0; i < pak_len; i++) {
        __m64 pak_val = _mm_set1_pi16(pak[i]);
        #pragma GCC unroll 4
        for (j = 0; j < 4; j++) {

            // out_pak[i*4 + j] = _mm_xor_si64(out_pak[i*4 + j], _mm_and_si64(_a[j], pak_val));

            __m64 res = _mm_and_si64(_a[j], pak_val);

            asm inline(
                "pxor %0, %1;"
                "movq %1, %0;"
            : "+m"(out_pak[i*4 + j]), "+y"(res)
            : 
            : );
        }
    }
}
#endif

#if (defined(__x86_64__) || defined(__i386__)) && defined(__SSE__)
void PERF_DEBUG_ATTRS __attribute__((visibility("default"))) fec_tx_col_op_sse(__m128* restrict out_pak, const unaligend_fec_int_t* restrict pak, size_t pak_len, fec_int_t a) {
    fec_idx_t i;

    int j;
    union {
        __m64 mmx[4];
        __m128 sse[2];
    } _a_tmp;
    const __m64 __a = _mm_set1_pi16(a);
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        __m64 res = _mm_srai_pi16(_mm_mullo_pi16(__a, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}), 15);
        asm inline(
                "movq %1, %0;"
            : "=m"(_a_tmp.mmx[j])
            : "y"(res)
            : );
    }

     __m128 _a_low;
     __m128 _a_high;

    asm inline(
                "movaps %1, %0;"
            : "=x"(_a_low)
            : "m"(_a_tmp.sse[0])
            : );
    
    asm inline(
                "movaps %1, %0;"
            : "=x"(_a_high)
            : "m"(_a_tmp.sse[1])
            : );

    for (i = 0; i < pak_len; i++) {

        uint16_t pak_val = pak[i];

        uint32_t pak_val_2 =  (pak_val | (((uint32_t)pak_val) << 16));
        __m128 pak_val_8;

        __m128 res;

        asm inline(
                "movss %[pak_val_2], %[pak_val_8];"
                "shufps $0, %[pak_val_8], %[pak_val_8];"
            : [pak_val_8] "=x"(pak_val_8)
            : [pak_val_2] "m"(pak_val_2)
            : );

        res = pak_val_8;

        asm inline(
                "andps %[_a_low], %[res];"
            : [res] "+x"(res)
            : [_a_low] "x"(_a_low)
            : );

        asm inline(
                "xorps %0, %1;"
                "movaps %1, %0;"
            : "+m"(out_pak[2*i]), "+x"(res)
            : 
            : );

        res = pak_val_8;

        asm inline(
                "andps %[_a_high], %[res];"
            : [res] "+x"(res)
            : [_a_high] "x"(_a_high)
            : );

        asm inline(
                "xorps %0, %1;"
                "movaps %1, %0;"
            : "+m"(out_pak[2*i + 1]), "+x"(res)
            : 
            : );

        // out_pak[2*i] = _mm_xor_ps(out_pak[2*i], _mm_and_ps(_a_low, pak_val_8));
        // out_pak[2*i + 1] = _mm_xor_ps(out_pak[2*i + 1], _mm_and_ps(_a_high, pak_val_8));
    }
}
#endif

void fec_tx_col_perf_to_norm_mmx(unaligend_fec_int_t* restrict out_pak, const uint16_t* restrict perf_pak, size_t pak_len) {
    fec_idx_t i;
    int j;
    for (i = 0; i < pak_len; i++) {

        uint32_t res = 0;
        for (j = 0; j < 16; j++) {
            res ^= ((uint32_t)perf_pak[i*16 + j]) << j;
        }

        uint16_t res4 = res;
        uint32_t carry = res >> 16;

        carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

        res4 ^= carry;

        carry = carry >> 16;
        
        res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

        out_pak[i] = res4;
    }
}

#endif

#if defined(FEC_HAS_64BIT)

void fec_tx_init_perf_arr_reg64(uint64_t* restrict out_pak, size_t pak_len) {
    memset(out_pak, 0, pak_len*sizeof(out_pak[0])*4);
}

void PERF_DEBUG_ATTRS fec_tx_col_op_reg64(uint64_t* restrict out_pak, const unaligend_fec_int_t* restrict pak, size_t pak_len, fec_int_t a) {

    int j;
    uint64_t _a[4];

    fec_idx_t i;
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        _a[j] = a * ((1ULL<<((16+1)*(4*j+0) - 64*j)) + (1ULL<<((16+1)*(4*j+1) - 64*j)) + (1ULL<<((16+1)*(4*j+2) - 64*j)) + (1ULL<<((16+1)*(4*j+3) - 64*j)));
        _a[j] >>= (16-1);
        _a[j] &= (1ULL<<(16*0)) | (1ULL<<(16*1)) | (1ULL<<(16*2)) | (1ULL<<(16*3));
        _a[j] *= (1<<16) - 1;
    }

    for (i = 0; i < pak_len; i++) {
        uint64_t b = (pak[i] * ((1ULL<<(16*0)) + (1ULL<<(16*1)) + (1ULL<<(16*2)) + (1ULL<<(16*3))));

        #pragma GCC unroll 4
        for (j = 0; j < 4; j++) {
            out_pak[i*4+j] ^= _a[j] & b;
        }
    }
}

void fec_tx_col_perf_to_norm_reg64(unaligend_fec_int_t* restrict out_pak, const uint16_t* restrict perf_pak, size_t pak_len) {
    fec_idx_t i;
    int j;
    for (i = 0; i < pak_len; i++) {

        uint32_t res = 0;
        for (j = 0; j < 16; j++) {
            res ^= ((uint32_t)perf_pak[i*16 + j]) << (16 - 1 - j);
        }

        uint16_t res4 = res;
        uint32_t carry = res >> 16;

        carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

        res4 ^= carry;

        carry = carry >> 16;
        
        res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

        out_pak[i] = res4;
    }
}

#endif


#ifdef FEC_HAS_32BIT

void fec_tx_init_perf_arr_reg32(uint32_t* restrict out_pak, size_t pak_len) {
    memset(out_pak, 0, pak_len*sizeof(out_pak[0])*8);
}

void PERF_DEBUG_ATTRS fec_tx_col_op_reg32(uint32_t* restrict out_pak, const unaligend_fec_int_t* restrict pak, size_t pak_len, fec_int_t a) {

    int j;
    uint32_t _a[8];

    fec_idx_t i;
    
    #pragma GCC unroll 8
    for (j = 0; j < 8; j++) {
        _a[j] = ((uint16_t)(((int16_t)(a << (16 - 1 - (j*2 + 0)))) >> 15)) | (((uint32_t)((uint16_t)(((int16_t)(a << (16 - 1 - (j*2 + 1)))) >> 15))) << 16);
    }

    for (i = 0; i < pak_len; i++) {
        uint32_t b = pak[i] | (((uint32_t)pak[i]) << 16);

        #pragma GCC unroll 8
        for (j = 0; j < 8; j++) {
            out_pak[i*8+j] ^= _a[j] & b;
        }
    }
}

void fec_tx_col_perf_to_norm_reg32(unaligend_fec_int_t* restrict out_pak, const uint16_t* restrict perf_pak, size_t pak_len) {
    fec_idx_t i;
    int j;
    for (i = 0; i < pak_len; i++) {

        uint32_t res = 0;
        for (j = 0; j < 16; j++) {
            res ^= ((uint32_t)perf_pak[i*16 + j]) << j;
        }

        uint16_t res4 = res;
        uint32_t carry = res >> 16;

        carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

        res4 ^= carry;

        carry = carry >> 16;
        
        res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

        out_pak[i] = res4;
    }
}

#endif


// TODO: idx can be fec_int_t
fec_status_t fec_tx_get_redundancy_pak(const fec_tx_state_t *tx_state, const fec_inv_cache_t *inv_cache, fec_idx_t idx, void *pak) {
    fec_idx_t n = tx_state->n;
    size_t pak_len = tx_state->pak_len;
    fec_idx_t i;
    size_t j;
    unaligend_fec_int_t* out_pak = (unaligend_fec_int_t*)pak;
    const unaligend_fec_int_t* restrict * restrict paks = tx_state->paks;
    
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_TX_OPT)
    fec_perf_int_t* restrict tmp_pak = tx_state->tmp_pak;
    tmp_pak = (fec_perf_int_t*)(((uintptr_t)tmp_pak + __alignof__(fec_perf_int_t) - 1) & (-__alignof__(fec_perf_int_t)));
    const fec_int_t* restrict inv_arr = inv_cache->inv_arr - 1;
#endif

    if (n == 0) {
        return FEC_STATUS_SUCCESS;
    }

    // subtract can't underflow cause n <= (((fec_idx_t)1)<<(sizeof(fec_int_t)*8))
    // idx = "k" - 1
    // so we check n + k - 1 <= MAX_FEC_INT + 1
    if (idx > (((fec_idx_t)1)<<(sizeof(fec_int_t)*8)) - n) {
        return FEC_STATUS_INVALID_PARAMS;
    }

    // inv_arr_sz is a power of 2 so this comparision is ok
    if (idx != 0 && inv_cache->inv_arr_sz < n + idx) {
        return FEC_STATUS_INV_CACHE_TOO_SMALL;
    }

    // TODO: this does not affect the time complexity of the function, but maybe remember number in struct
    for (i = 0; i < n; i++) {
        if (tx_state->paks[i] == NULL) {
            return FEC_STATUS_MORE_PACKETS_NEEDED;
        }
    }

    // for (j = 0; j < pak_len; j++) {
    //     fec_int_t res = 0;
    //     for (i = 0; i < n; i++) {
    //         if(idx == 0) {
    //             res = poly_add(res, tx_state->paks[i][j]);
    //         } else {
    //             fec_int_t a_i = poly_add(n + idx - 1, i);
    //             res = poly_add(res, poly_mul(tx_state->paks[i][j], _fec_inv(inv_cache, a_i)));
    //         }
    //     }
    //     out_pak[j] = res;
    // }

    memset(out_pak, 0, pak_len*sizeof(fec_int_t));

    if (idx == 0) {
        for (i = 0; i < n; i++) {
            const unaligend_fec_int_t* pak = paks[i];
            for (j = 0; j < pak_len; j++) {
                out_pak[j] = poly_add(out_pak[j], pak[j]);
            }
        }
        return FEC_STATUS_SUCCESS;
    }

    // for (j = 0; j < pak_len; j++) {
    //     out_pak[j] = rand();
    // }

#if defined(FEC_HAS_CLMUL64)

    // seems to be bound by L3 cache size for performance

    fec_tx_init_perf_arr(tmp_pak, pak_len);
    for (i = 0; i < n; i++) {
        fec_int_t a_i = inv_arr[poly_add(n + idx - 1, i)];
        const unaligend_fec_int_t* pak = paks[i];
        fec_tx_col_op(tmp_pak, pak, pak_len, a_i);
    }
    fec_tx_col_perf_to_norm(out_pak, tmp_pak, pak_len);
    // size_t k;
    
    // for (k = 0; k < pak_len; k += 512*2048U) {
    //     size_t cur_pak_len = MIN(pak_len - k, 512*2048U);
    //     fec_tx_init_perf_arr(tmp_pak, cur_pak_len);
    //     for (i = 0; i < n; i++) {
    //         fec_int_t a_i = inv_arr[poly_add(n + idx - 1, i)];
    //         const unaligend_fec_int_t* pak = paks[i];
    //         // for (j = 0; j < pak_len; j++) {
    //         //     // if (idx == 0) {
    //         //     //     out_pak[j] = poly_add(out_pak[j], pak[j]);
    //         //     // } else {
    //         //         out_pak[j] = poly_add(out_pak[j], poly_mul(pak[j], a_i));
    //         //     // }
    //         // }
    //         fec_tx_col_op(tmp_pak, &pak[k], cur_pak_len, a_i);
    //     }
    //     fec_tx_col_perf_to_norm(&out_pak[k], tmp_pak, cur_pak_len);
    // }
#elif defined(FEC_HAS_128_INT_VEC) || defined(FEC_HAS_64_INT_VEC)
    size_t k;
    fec_idx_t m;

    fec_tx_init_perf_arr_avx(tmp_pak, pak_len);
    for (i = 0; i < n; i++) {
        fec_int_t a_i = inv_arr[poly_add(n + idx - 1, i)];
        const unaligend_fec_int_t* pak = paks[i];
#if defined(FEC_HAS_128_INT_VEC) || (defined(FEC_HAS_64_INT_VEC) && !(defined(__x86_64__) || defined(__i386__)))
        fec_tx_col_op_avx(tmp_pak, pak, pak_len, a_i);
#elif (defined(__x86_64__) || defined(__i386__)) && defined(__SSE__)
        fec_tx_col_op_sse((__m128*)tmp_pak, pak, pak_len, a_i);
#elif (defined(__x86_64__) || defined(__i386__)) && defined(__MMX__)
        fec_tx_col_op_mmx((__m64*)tmp_pak, pak, pak_len, a_i);
#endif
    }

#if defined(FEC_HAS_128_INT_VEC)
    fec_tx_col_perf_to_norm_avx(out_pak, tmp_pak, pak_len);
#else
    fec_tx_col_perf_to_norm_mmx(out_pak, (const uint16_t*)tmp_pak, pak_len);
#endif
    
#if !defined(FEC_HAS_128_INT_VEC) && defined(FEC_HAS_64_INT_VEC) && (defined(__x86_64__) || defined(__i386__))
    _mm_empty();
#endif

    // // seems to be bound by L2 cache size for performance
    
    // #define MULT 6

    // for (k = 0; k < pak_len; k += ((size_t)1024U*MULT)) {
    //     size_t cur_pak_len = MIN(pak_len - k, ((size_t)1024U*MULT));
    //     fec_tx_init_perf_arr_avx(tmp_pak, cur_pak_len);
    //     for (i = 0; i < n; i++) {
    //         fec_int_t a_i = inv_arr[poly_add(n + idx - 1, i)];
    //         const unaligend_fec_int_t* pak = paks[i];
    //         // for (j = 0; j < pak_len; j++) {
    //         //     // if (idx == 0) {
    //         //     //     out_pak[j] = poly_add(out_pak[j], pak[j]);
    //         //     // } else {
    //         //         out_pak[j] = poly_add(out_pak[j], poly_mul(pak[j], a_i));
    //         //     // }
    //         // }
    //         fec_tx_col_op_avx(tmp_pak, &pak[k], cur_pak_len, a_i);
    //     }
    //     fec_tx_col_perf_to_norm_avx(&out_pak[k], tmp_pak, cur_pak_len);
    // }

// #ifndef N_BLOCK
//     #define N_BLOCK 1
// #endif
//     #define LEN_BLOCK (32*1024 - 4*1024 - 8*N_BLOCK) / (64 + 2*N_BLOCK)
//     //#define LEN_BLOCK 1024

//     fec_tx_init_perf_arr_avx(tmp_pak, pak_len);

//     for(m = 0; m < n; m += N_BLOCK) {
//         fec_idx_t cur_n = MIN(n - m, N_BLOCK);

//         for (k = 0; k < pak_len; k += LEN_BLOCK) {

//             // size_t ii;
//             // for (ii = 0; ii < LEN_BLOCK * sizeof(tmp_pak[0]); ii += 0x40) {
//             //     _mm_prefetch(((char*)&tmp_pak[k + LEN_BLOCK]) + ii, _MM_HINT_ET0);
//             // }
            
//             size_t cur_pak_len = MIN(pak_len - k, LEN_BLOCK);
            
//             for (i = 0; i < cur_n; i++) {
//                 fec_int_t a_i = inv_arr[poly_add(n + idx - 1, m + i)];
//                 const unaligend_fec_int_t* pak = paks[m + i];
//                 // for (j = 0; j < pak_len; j++) {
//                 //     // if (idx == 0) {
//                 //     //     out_pak[j] = poly_add(out_pak[j], pak[j]);
//                 //     // } else {
//                 //         out_pak[j] = poly_add(out_pak[j], poly_mul(pak[j], a_i));
//                 //     // }
//                 // }
//                 fec_tx_col_op_avx(&tmp_pak[k], &pak[k], cur_pak_len, a_i);
//             }
            
//         }
//     }

//     fec_tx_col_perf_to_norm_avx(out_pak, tmp_pak, pak_len);

#elif defined(FEC_HAS_64BIT)
    // fec_tx_init_perf_arr_reg64((uint64_t*)tmp_pak, pak_len);
    // for (i = 0; i < n; i++) {
    //     fec_int_t a_i = _fec_inv(inv_cache, poly_add(n + idx - 1, i));
    //     const unaligend_fec_int_t* pak = paks[i];
    //     fec_tx_col_op_reg64((uint64_t*)tmp_pak, pak, pak_len, a_i);
    // }
    // fec_tx_col_perf_to_norm_reg64(out_pak, (const uint16_t*)tmp_pak, pak_len);

    size_t k;

    // 64-90 is also sweet spot
    // but seems like L3 bound
    
    for (k = 0; k < pak_len; k += 8*1024U) {
        size_t cur_pak_len = MIN(pak_len - k, 8*1024U);
        fec_tx_init_perf_arr_reg64((uint64_t*)tmp_pak, cur_pak_len);
        for (i = 0; i < n; i++) {
            fec_int_t a_i = inv_arr[poly_add(n + idx - 1, i)];
            const unaligend_fec_int_t* pak = paks[i];
            // for (j = 0; j < pak_len; j++) {
            //     // if (idx == 0) {
            //     //     out_pak[j] = poly_add(out_pak[j], pak[j]);
            //     // } else {
            //         out_pak[j] = poly_add(out_pak[j], poly_mul(pak[j], a_i));
            //     // }
            // }
            fec_tx_col_op_reg64((uint64_t*)tmp_pak, &pak[k], cur_pak_len, a_i);
        }
        fec_tx_col_perf_to_norm_reg64(&out_pak[k], (const uint16_t*)tmp_pak, cur_pak_len);
    }
#elif defined(FEC_HAS_32BIT)
    // fec_tx_init_perf_arr_reg32((uint32_t*)tmp_pak, pak_len);
    // for (i = 0; i < n; i++) {
    //     fec_int_t a_i = _fec_inv(inv_cache, poly_add(n + idx - 1, i));
    //     const unaligend_fec_int_t* pak = paks[i];
    //     fec_tx_col_op_reg32((uint32_t*)tmp_pak, pak, pak_len, a_i);
    // }
    // fec_tx_col_perf_to_norm_reg32(out_pak, (const uint16_t*)tmp_pak, pak_len);

    size_t k;
    
    for (k = 0; k < pak_len; k += 8*1024U) {
        size_t cur_pak_len = MIN(pak_len - k, 8*1024U);
        fec_tx_init_perf_arr_reg32((uint32_t*)tmp_pak, cur_pak_len);
        for (i = 0; i < n; i++) {
            fec_int_t a_i = inv_arr[poly_add(n + idx - 1, i)];
            const unaligend_fec_int_t* pak = paks[i];
            fec_tx_col_op_reg32((uint32_t*)tmp_pak, &pak[k], cur_pak_len, a_i);
        }
        fec_tx_col_perf_to_norm_reg32(&out_pak[k], (const uint16_t*)tmp_pak, cur_pak_len);
    }
#else
    for (i = 0; i < n; i++) {
        fec_int_t a_i = _fec_inv(inv_cache, poly_add(n + idx - 1, i));
        const unaligend_fec_int_t* pak = paks[i];
        for (j = 0; j < pak_len; j++) {
            // if (idx == 0) {
            //     out_pak[j] = poly_add(out_pak[j], pak[j]);
            // } else {
                out_pak[j] = poly_add(out_pak[j], poly_mul(pak[j], a_i));
            // }
        }
    }
#endif

    return FEC_STATUS_SUCCESS;
}

#ifndef FEC_MIN_MEM
#ifndef FEC_LARGE_K
fec_status_t fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state, const fec_inv_cache_t *inv_cache) {
    fec_idx_t n = rx_state->n;
    size_t pak_len = rx_state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present = 0;

    fec_idx_t i, j;
    fec_idx_t y_i, y_j;
    size_t ii;

    if (rx_state->num_redundant < num_y_missing) {
        return FEC_STATUS_MORE_PACKETS_NEEDED;
    }

    if (num_y_missing == 0) {
        return FEC_STATUS_SUCCESS;
    }

    // n == 0 -> num_y_missing == 0
    // k == 0 -> no redundants and can recover -> num_y_missing == 0
    // here n,k > 0
    // max_x = n + k - 2 -> n + k - 1 = max_x + 1
    // inv_cache->inv_arr_sz < n + k - 1 <-> inv_cache->inv_arr_sz <= max_x
    if (inv_cache->inv_arr_sz <= rx_state->max_x) {
        return FEC_STATUS_INV_CACHE_TOO_SMALL;
    }

    has_one_row = rx_state->redundancy_paks[0] != NULL;
    if (num_y_missing >= has_one_row) {
        num_x_present = num_y_missing - has_one_row;
    }

    for (y_i = 0, i = 0, j = 0; y_i < n; y_i++) {
        if (rx_state->info_paks[y_i] != NULL) {
            rx_state->present_y[j] = y_i;
            j++;
        } else {
            rx_state->missing_y[i] = y_i;
            i++;
        }
    }

    for (i = 0; i < num_x_present; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_y_missing; j++) {
            res = poly_mul(res, poly_add(rx_state->present_x[i], rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            res = poly_mul(res, _fec_inv(inv_cache, poly_add(rx_state->present_x[i], rx_state->present_x[j])));
        }
        rx_state->pi_xy_div_xx[i] = res;
    }

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_x_present; j++) {
            res = poly_mul(res, poly_add(rx_state->missing_y[i], rx_state->present_x[j]));
        }
        for (j = 0; j < num_y_missing; j++) {
            if(j == i) {
                continue;
            }
            res = poly_mul(res, _fec_inv(inv_cache, poly_add(rx_state->missing_y[i], rx_state->missing_y[j])));
        }
        rx_state->pi_yx_div_yy[i] = res;
    }

    for (y_i = 0, i = 0; y_i < n; y_i++) {
        if (rx_state->info_paks[y_i] == NULL) {
            continue;
        }
        fec_int_t res = 1;
        for (j = 0; j < num_y_missing; j++) {
            res = poly_mul(res, poly_add(y_i, rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            res = poly_mul(res, _fec_inv(inv_cache, poly_add(y_i, rx_state->present_x[j])));
        }
        rx_state->pi_ycomp_y_div_ycomp_x[i] = res;
        i++;
    }

    for (ii = 0; ii < pak_len; ii++) {
        fec_int_t res;
        fec_int_t tmp_vec_1s = 0;
        for (y_i = 0, i = 0; y_i < n; y_i++) {
            if (rx_state->info_paks[y_i] == NULL) {
                continue;
            }
            res = poly_mul(rx_state->info_paks[y_i][ii], rx_state->pi_ycomp_y_div_ycomp_x[i]);
            rx_state->tmp_vec_info[i] = res;
            i++;
        }
        if(has_one_row) {
            tmp_vec_1s = rx_state->redundancy_paks[0][ii];
        }
        //for (x_i = 0, i = 0; x_i < k - 1; x_i++) {
        for(i = 0; i < num_x_present; i++) {
            res = poly_mul(rx_state->redundancy_paks[rx_state->present_x[i] - n + 1][ii], rx_state->pi_xy_div_xx[i]);
            rx_state->tmp_vec_redundancy[i] = res;
        }

        for (i = 0; i < num_y_missing; i++) {

            res = 0;

            for (y_j = 0, j = 0; y_j < n; y_j++) {
                if (rx_state->info_paks[y_j] == NULL) {
                    continue;
                }
                res = poly_add(res, poly_mul(rx_state->tmp_vec_info[j], _fec_inv(inv_cache, poly_add(y_j, rx_state->missing_y[i]))));
                j++;
            }
            if (has_one_row) {
                res = poly_add(res, tmp_vec_1s);
            }
            //for (x_j = 0, j = 0; x_j < k - 1; x_j++) {
            for(j = 0; j < num_x_present; j++) {
                res = poly_add(res, poly_mul(rx_state->tmp_vec_redundancy[j], _fec_inv(inv_cache, poly_add(rx_state->present_x[j], rx_state->missing_y[i]))));
            }

            res = poly_mul(res, rx_state->pi_yx_div_yy[i]);

            if (i == 0 && has_one_row) {
                rx_state->redundancy_paks[0][ii] = res;
            } else {
                rx_state->redundancy_paks[rx_state->present_x[i - has_one_row] - n + 1][ii] = res;
            }
        }

    }

    if (has_one_row) {
        rx_state->info_paks[rx_state->missing_y[0]] = rx_state->redundancy_paks[0];
    }
    for (i = 0; i < num_x_present; i++) {
        rx_state->info_paks[rx_state->missing_y[i + has_one_row]] = rx_state->redundancy_paks[rx_state->present_x[i] - n + 1];
    }

    return FEC_STATUS_SUCCESS;
}

#else

fec_status_t fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state, const fec_inv_cache_t *inv_cache) {
    fec_idx_t n = rx_state->n;
    size_t pak_len = rx_state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present = 0;
    fec_int_t *present_x;

    fec_idx_t i, j;
    fec_idx_t y_i, y_j;
    size_t ii;

    if (rx_state->num_info + rx_state->num_redundant + (rx_state->ones_pak != NULL) < n) {
        return FEC_STATUS_MORE_PACKETS_NEEDED;
    }

    if (num_y_missing == 0) {
        return FEC_STATUS_SUCCESS;
    }

    // n == 0 -> num_y_missing == 0
    // k == 0 -> no redundants and can recover -> num_y_missing == 0
    // here n,k > 0
    // max_x = n + k - 2 -> n + k - 1 = max_x + 1
    // inv_cache->inv_arr_sz < n + k - 1 <-> inv_cache->inv_arr_sz <= max_x
    if (inv_cache->inv_arr_sz <= rx_state->max_x) {
        return FEC_STATUS_INV_CACHE_TOO_SMALL;
    }

    has_one_row = (rx_state->ones_pak != NULL);
    num_x_present = num_y_missing - has_one_row;

    present_x = &rx_state->pak_xy_arr[n - num_x_present];

    for (y_i = 0, i = 0; i < num_y_missing; y_i++) {
        if ((rx_state->received_paks_bitmap[y_i / 8] & (1<<(y_i & (8-1)))) == 0) {
            rx_state->missing_y[i] = y_i;
            i++;
        }
    }

    for (i = 0; i < num_x_present; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_y_missing; j++) {
            res = poly_mul(res, poly_add(present_x[i], rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            res = poly_mul(res, _fec_inv(inv_cache, poly_add(present_x[i], present_x[j])));
        }
        rx_state->pi_xy_div_xx[i] = res;
    }

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_x_present; j++) {
            res = poly_mul(res, poly_add(rx_state->missing_y[i], present_x[j]));
        }
        for (j = 0; j < num_y_missing; j++) {
            if(j == i) {
                continue;
            }
            res = poly_mul(res, _fec_inv(inv_cache, poly_add(rx_state->missing_y[i], rx_state->missing_y[j])));
        }
        rx_state->pi_yx_div_yy[i] = res;
    }

    for(i = 0; i < rx_state->num_info; i++) {
        y_i = rx_state->pak_xy_arr[i];
        fec_int_t res = 1;
        for (j = 0; j < num_y_missing; j++) {
            res = poly_mul(res, poly_add(y_i, rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            res = poly_mul(res, _fec_inv(inv_cache, poly_add(y_i, present_x[j])));
        }
        rx_state->pi_ycomp_y_div_ycomp_x[i] = res;
    }

    for (ii = 0; ii < pak_len; ii++) {
        fec_int_t res;
        fec_int_t tmp_vec_1s = 0;
        for(i = 0; i < rx_state->num_info; i++) {
            y_i = rx_state->pak_xy_arr[i];
            res = poly_mul(rx_state->pak_arr[i][ii], rx_state->pi_ycomp_y_div_ycomp_x[i]);
            rx_state->tmp_vec_info[i] = res;
        }
        if(has_one_row) {
            tmp_vec_1s = rx_state->ones_pak[ii];
        }
        for(i = 0; i < num_x_present; i++) {
            res = poly_mul(rx_state->pak_arr[n - num_x_present + i][ii], rx_state->pi_xy_div_xx[i]);
            rx_state->tmp_vec_redundancy[i] = res;
        }

        for (i = 0; i < num_y_missing; i++) {

            res = 0;

            for(j = 0; j < rx_state->num_info; j++) {
                y_j = rx_state->pak_xy_arr[j];
                res = poly_add(res, poly_mul(rx_state->tmp_vec_info[j], _fec_inv(inv_cache, poly_add(y_j, rx_state->missing_y[i]))));
            }
            if (has_one_row) {
                res = poly_add(res, tmp_vec_1s);
            }
            for(j = 0; j < num_x_present; j++) {
                res = poly_add(res, poly_mul(rx_state->tmp_vec_redundancy[j], _fec_inv(inv_cache, poly_add(present_x[j], rx_state->missing_y[i]))));
            }

            res = poly_mul(res, rx_state->pi_yx_div_yy[i]);

            if (i == 0 && has_one_row) {
                rx_state->ones_pak[ii] = res;
            } else {
                rx_state->pak_arr[n - num_x_present + i - has_one_row][ii] = res;
            }
        }

    }

    if (has_one_row) {
        rx_state->pak_xy_arr[n - num_x_present - 1] = rx_state->missing_y[0];
        rx_state->pak_arr[n - num_x_present - 1] = rx_state->ones_pak;
    }
    for (i = 0; i < num_x_present; i++) {
        rx_state->pak_xy_arr[n - num_x_present + i] = rx_state->missing_y[i + has_one_row];
    }

    for (i = 0; i < n; i++) {
        while (rx_state->pak_xy_arr[i] != i) {
            fec_int_t idx = rx_state->pak_xy_arr[i];

            {
                fec_int_t tmp;
                tmp = rx_state->pak_xy_arr[idx];
                rx_state->pak_xy_arr[idx] = idx;
                rx_state->pak_xy_arr[i] = tmp;
            }

            {
                unaligend_fec_int_t *tmp;
                tmp = rx_state->pak_arr[idx];
                rx_state->pak_arr[idx] = rx_state->pak_arr[i];
                rx_state->pak_arr[i] = tmp;
            }
        }
    }

    return FEC_STATUS_SUCCESS;
}

#endif

#else

#ifndef FEC_LARGE_K

fec_status_t fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state, const fec_inv_cache_t *inv_cache) {
    fec_idx_t n = rx_state->n;
    size_t pak_len = rx_state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present = 0;

    fec_idx_t i, j;
    fec_idx_t y_i, y_j;
    size_t ii;

    if (rx_state->num_redundant < num_y_missing) {
        return FEC_STATUS_MORE_PACKETS_NEEDED;
    }

    if (num_y_missing == 0) {
        return FEC_STATUS_SUCCESS;
    }

    // n == 0 -> num_y_missing == 0
    // k == 0 -> no redundants and can recover -> num_y_missing == 0
    // here n,k > 0
    // max_x = n + k - 2 -> n + k - 1 = max_x + 1
    // inv_cache->inv_arr_sz < n + k - 1 <-> inv_cache->inv_arr_sz <= max_x
    if (inv_cache->inv_arr_sz <= rx_state->max_x) {
        return FEC_STATUS_INV_CACHE_TOO_SMALL;
    }

    has_one_row = rx_state->redundancy_paks[0] != NULL;
    if (num_y_missing >= has_one_row) {
        num_x_present = num_y_missing - has_one_row;
    }

    for (y_i = 0, i = 0; i < num_y_missing; y_i++) {
        if (rx_state->info_paks[y_i] == NULL) {
            rx_state->missing_y[i] = y_i;
            i++;
        }
    }

    for (y_i = 0, i = 0; y_i < n; y_i++) {
        if (rx_state->info_paks[y_i] == NULL) {
            continue;
        }
        fec_int_t pi_ycomp_y_div_ycomp_x_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, poly_add(y_i, rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, _fec_inv(inv_cache, poly_add(y_i, rx_state->present_x[j])));
        }
        for (ii = 0; ii < pak_len; ii++) {
            rx_state->info_paks[y_i][ii] = poly_mul(rx_state->info_paks[y_i][ii], pi_ycomp_y_div_ycomp_x_i);
        }
        i++;
    }

    for (i = 0; i < num_x_present; i++) {
        fec_int_t pi_xy_div_xx_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, poly_add(rx_state->present_x[i], rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, _fec_inv(inv_cache, poly_add(rx_state->present_x[i], rx_state->present_x[j])));
        }

        unaligend_fec_int_t *pak = rx_state->redundancy_paks[rx_state->present_x[i] - n + 1];
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_xy_div_xx_i);
        }
    }

    for (ii = 0; ii < pak_len; ii++) {
        fec_int_t res;

        for (i = 0; i < num_y_missing; i++) {

            res = 0;

            for (y_j = 0, j = 0; y_j < n; y_j++) {
                if (rx_state->info_paks[y_j] == NULL) {
                    continue;
                }
                res = poly_add(res, poly_mul(rx_state->info_paks[y_j][ii], _fec_inv(inv_cache, poly_add(y_j, rx_state->missing_y[i]))));
                j++;
            }
            if (has_one_row) {
                res = poly_add(res, rx_state->redundancy_paks[0][ii]);
            }
            //for (x_j = 0, j = 0; x_j < k - 1; x_j++) {
            for(j = 0; j < num_x_present; j++) {
                res = poly_add(res, poly_mul(rx_state->redundancy_paks[rx_state->present_x[j] - n + 1][ii], _fec_inv(inv_cache, poly_add(rx_state->present_x[j], rx_state->missing_y[i]))));
            }

            rx_state->tmp_recovered_ints[i] = res;
        }

        for (i = 0; i < num_y_missing; i++) {
            if (i == 0 && has_one_row) {
                rx_state->redundancy_paks[0][ii] = rx_state->tmp_recovered_ints[i];
            } else {
                rx_state->redundancy_paks[rx_state->present_x[i - has_one_row] - n + 1][ii] = rx_state->tmp_recovered_ints[i];
            }
        }
    }

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t pi_yx_div_yy_i = 1;
        for (j = 0; j < num_x_present; j++) {
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, poly_add(rx_state->missing_y[i], rx_state->present_x[j]));
        }
        for (j = 0; j < num_y_missing; j++) {
            if(j == i) {
                continue;
            }
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, _fec_inv(inv_cache, poly_add(rx_state->missing_y[i], rx_state->missing_y[j])));
        }

        for (ii = 0; ii < pak_len; ii++) {
            if (i == 0 && has_one_row) {
                rx_state->redundancy_paks[0][ii] = poly_mul(rx_state->redundancy_paks[0][ii], pi_yx_div_yy_i);
            } else {
                rx_state->redundancy_paks[rx_state->present_x[i - has_one_row] - n + 1][ii] = poly_mul(rx_state->redundancy_paks[rx_state->present_x[i - has_one_row] - n + 1][ii], pi_yx_div_yy_i);
            }
        }
    }

    for (y_i = 0, i = 0; y_i < n; y_i++) {
        if (rx_state->info_paks[y_i] == NULL) {
            continue;
        }
        fec_int_t inv_pi_ycomp_y_div_ycomp_x_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i, _fec_inv(inv_cache, poly_add(y_i, rx_state->missing_y[j])));
        }
        for (j = 0; j < num_x_present; j++) {
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i,  poly_add(y_i, rx_state->present_x[j]));
        }
        for (ii = 0; ii < pak_len; ii++) {
            rx_state->info_paks[y_i][ii] = poly_mul(rx_state->info_paks[y_i][ii], inv_pi_ycomp_y_div_ycomp_x_i);
        }
        i++;
    }

    if (has_one_row) {
        rx_state->info_paks[rx_state->missing_y[0]] = rx_state->redundancy_paks[0];
    }
    for (i = 0; i < num_x_present; i++) {
        rx_state->info_paks[rx_state->missing_y[i + has_one_row]] = rx_state->redundancy_paks[rx_state->present_x[i] - n + 1];
    }

    return FEC_STATUS_SUCCESS;
}

#else

#include <emmintrin.h>

static void PERF_DEBUG_ATTRS __fec_rx_col_op(fec_perf_int_t* restrict recovered, const fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_int_t* restrict inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {

    inv_arr = inv_arr - 1;

    /*{
        fec_idx_t i;

        u32x4 _b = {pak_val, 0, 0, 0};
        
        for (i = 0; i < num_y_missing - 1; i+=2) {
            // uint32_t h = ((uint32_t*)missing_y)[i/2];
            // uint16_t a1 = inv_arr[poly_add(pak_xy, h)];
            // uint16_t a2 = inv_arr[poly_add(pak_xy, h >> 16)];
            uint16_t a1 = inv_arr[poly_add(pak_xy, missing_y[i])];
            uint16_t a2 = inv_arr[poly_add(pak_xy, missing_y[i+1])];

            u32x4 _a = {a1, a2, 0, 0};

            u32x4 _c = (u32x4)_mm_clmulepi64_si128((__m128i)_a, (__m128i)_b, 0);
            
            
            recovered[i] ^= _c[0];
            recovered[i+1] ^= _c[1];
        }

        if (i == num_y_missing - 1) {
            uint16_t a1 = inv_arr[poly_add(pak_xy, missing_y[i])];
            u32x4 _a = {a1, 0, 0, 0};
            u32x4 _c = (u32x4)_mm_clmulepi64_si128((__m128i)_a, (__m128i)_b, 0);
            recovered[i] ^= _c[0];
        }
    }*/

    // {
    //     fec_idx_t i;

    //     u16x16 shifts = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    //     u16x16 _b = (u16x16)(((i16x16)(pak_val << shifts)) >> 15);

    //     for (i = 0; i < num_y_missing; i++) {
    //         uint16_t a = inv_arr[poly_add(pak_xy, missing_y[i])];
    //         u16x16 _c = a & _b;
    //         recovered[i] ^= _c;
    //     }
    // }

    return;
/*
#if defined(_FEC_USE_POLY_MUL)
    fec_idx_t i;
    for (i = 0; i < num_y_missing; i++) {
        recovered[i] ^= poly_mul(pak_val, inv_arr[poly_add(pak_xy, missing_y[i])]);
    }
#elif defined(_FEC_USE_POLY_MUL16)
    fec_idx_t i;
    u16x16 aa = {pak_val,pak_val,pak_val,pak_val,pak_val,pak_val,pak_val,pak_val,
                pak_val,pak_val,pak_val,pak_val,pak_val,pak_val,pak_val,pak_val};
    for (i = 0; i < num_y_missing; i+=16) {
        u16x16 bb;

        // fec_idx_t j;
        // fec_idx_t num_to_copy = MIN(16, num_y_missing - i);
        // for(j = 0; j < num_to_copy; j++) {
        //     bb[j] = inv_arr[poly_add(pak_xy, missing_y[i+j])];
        // }

        uint16_t __attribute__((aligned(32))) tmp_arr[16];
        fec_idx_t j;
        fec_idx_t num_to_copy = MIN(16U, num_y_missing - i);
        for(j = 0; j < num_to_copy; j++) {
            tmp_arr[j] = inv_arr[poly_add(pak_xy, missing_y[i+j])];
        }
        bb = (u16x16)_mm256_load_si256((__m256i*)tmp_arr);
        

        // bb[0] = inv_arr[poly_add(pak_xy, missing_y[i+0])];
        // bb[1] = inv_arr[poly_add(pak_xy, missing_y[i+1])];
        // bb[2] = inv_arr[poly_add(pak_xy, missing_y[i+2])];
        // bb[3] = inv_arr[poly_add(pak_xy, missing_y[i+3])];
        // bb[4] = inv_arr[poly_add(pak_xy, missing_y[i+4])];
        // bb[5] = inv_arr[poly_add(pak_xy, missing_y[i+5])];
        // bb[6] = inv_arr[poly_add(pak_xy, missing_y[i+6])];
        // bb[7] = inv_arr[poly_add(pak_xy, missing_y[i+7])];
        // bb[8] = inv_arr[poly_add(pak_xy, missing_y[i+8])];
        // bb[9] = inv_arr[poly_add(pak_xy, missing_y[i+9])];
        // bb[10] = inv_arr[poly_add(pak_xy, missing_y[i+10])];
        // bb[11] = inv_arr[poly_add(pak_xy, missing_y[i+11])];
        // bb[12] = inv_arr[poly_add(pak_xy, missing_y[i+12])];
        // bb[13] = inv_arr[poly_add(pak_xy, missing_y[i+13])];
        // bb[14] = inv_arr[poly_add(pak_xy, missing_y[i+14])];
        // bb[15] = inv_arr[poly_add(pak_xy, missing_y[i+15])];
        u16x16 cc = poly_mul16(bb, aa);
        //u16x16 cc = poly_mul16_single(bb, pak_val);

        _mm256_storeu_si256((__m256i_u*)&recovered[i],  (__m256i)((u16x16)_mm256_loadu_si256((__m256i_u*)&recovered[i]) ^ cc));

        // recovered[i+0] ^= cc[0];
        // recovered[i+1] ^= cc[1];
        // recovered[i+2] ^= cc[2];
        // recovered[i+3] ^= cc[3];
        // recovered[i+4] ^= cc[4];
        // recovered[i+5] ^= cc[5];
        // recovered[i+6] ^= cc[6];
        // recovered[i+7] ^= cc[7];
        // recovered[i+8] ^= cc[8];
        // recovered[i+9] ^= cc[9];
        // recovered[i+10] ^= cc[10];
        // recovered[i+11] ^= cc[11];
        // recovered[i+12] ^= cc[12];
        // recovered[i+13] ^= cc[13];
        // recovered[i+14] ^= cc[14];
        // recovered[i+15] ^= cc[15];
    }
#elif defined(_FEC_USE_POLY_MUL8)
    fec_idx_t i;
    u16x8 aa = {pak_val,pak_val,pak_val,pak_val,pak_val,pak_val,pak_val,pak_val};
    // union {
    //     uint16_t __attribute__((aligned(16))) tmp_arr[8];
    //     u16x8 _bb;
    // } u;
    
    for (i = 0; i < num_y_missing; i+=8) {
        u16x8 bb;

        // bb = (u16x8)_mm_loadu_si128((__m128i_u*)&missing_y[i]);
        // bb ^= pak_xy;
        // bb[0] = inv_arr[bb[0]];
        // bb[1] = inv_arr[bb[1]];
        // bb[2] = inv_arr[bb[2]];
        // bb[3] = inv_arr[bb[3]];
        // bb[4] = inv_arr[bb[4]];
        // bb[5] = inv_arr[bb[5]];
        // bb[6] = inv_arr[bb[6]];
        // bb[7] = inv_arr[bb[7]];

        // int j;
        // for(j = 0; j < 8; j++) {
        //     bb[j] = inv_arr[poly_add(pak_xy, missing_y[i+j])];
        // }

        // int j;
        // for(j = 0; j < 8; j++) {
        //     u.tmp_arr[j] = inv_arr[poly_add(pak_xy, missing_y[i+j])];
        // }
        // bb = u._bb;
        uint16_t __attribute__((aligned(16))) tmp_arr[8];
        fec_idx_t j;
        fec_idx_t num_to_copy = MIN(8U, num_y_missing - i);
        for(j = 0; j < num_to_copy; j++) {
            tmp_arr[j] = inv_arr[poly_add(pak_xy, missing_y[i+j])];
        }
        bb = (u16x8)_mm_load_si128((__m128i*)tmp_arr);

        // bb[0] = inv_arr[poly_add(pak_xy, missing_y[i+0])];
        // bb[1] = inv_arr[poly_add(pak_xy, missing_y[i+1])];
        // bb[2] = inv_arr[poly_add(pak_xy, missing_y[i+2])];
        // bb[3] = inv_arr[poly_add(pak_xy, missing_y[i+3])];
        // bb[4] = inv_arr[poly_add(pak_xy, missing_y[i+4])];
        // bb[5] = inv_arr[poly_add(pak_xy, missing_y[i+5])];
        // bb[6] = inv_arr[poly_add(pak_xy, missing_y[i+6])];
        // bb[7] = inv_arr[poly_add(pak_xy, missing_y[i+7])];
        u16x8 cc = poly_mul8(aa, bb);

        _mm_storeu_si128((__m128i_u*)&recovered[i],  (__m128i)((u16x8)_mm_loadu_si128((__m128i_u*)&recovered[i]) ^ cc));

        // recovered[i+0] ^= cc[0];
        // recovered[i+1] ^= cc[1];
        // recovered[i+2] ^= cc[2];
        // recovered[i+3] ^= cc[3];
        // recovered[i+4] ^= cc[4];
        // recovered[i+5] ^= cc[5];
        // recovered[i+6] ^= cc[6];
        // recovered[i+7] ^= cc[7];
    }
#elif defined(_FEC_USE_POLY_MUL4_MMX)
    fec_idx_t i;
    __m64 aa = (__m64)((u16x4){pak_val,pak_val,pak_val,pak_val});
    
    for (i = 0; i < num_y_missing; i+=4) {
        __m64 bb;

        uint16_t __attribute__((aligned(8))) tmp_arr[4];
        fec_idx_t j;
        fec_idx_t num_to_copy = MIN(4U, num_y_missing - i);
        for(j = 0; j < num_to_copy; j++) {
            tmp_arr[j] = inv_arr[poly_add(pak_xy, missing_y[i+j])];
        }
        bb = *((__m64*)tmp_arr);

        __m64 cc = poly_mul4_mmx(aa, bb);

        *((__m64*)&recovered[i]) = _mm_xor_si64(*((__m64*)&recovered[i]), cc);
    }
#elif defined(_FEC_USE_POLY_MUL4)
    fec_idx_t i;
    uint64_t aa = ((uint64_t)pak_val << 0) | ((uint64_t)pak_val << 16) | ((uint64_t)pak_val << 32) | ((uint64_t)pak_val << 48);
    for (i = 0; i < num_y_missing; i+=4) {
        uint64_t bb = 0;
        fec_idx_t j;
        fec_idx_t num_to_copy = MIN(4U, num_y_missing - i);
        for (j = 0; j < num_to_copy; j++) {
            bb |= ((uint64_t)inv_arr[poly_add(pak_xy, missing_y[i+j])]) << (j*16);
        }
        uint64_t cc = poly_mul4(aa, bb);

        *((uint64_t*)&recovered[i]) ^= cc;
        // recovered[i+0] ^= (uint16_t)(cc >> 0);
        // recovered[i+1] ^= (uint16_t)(cc >> 16);
        // recovered[i+2] ^= (uint16_t)(cc >> 32);
        // recovered[i+3] ^= (uint16_t)(cc >> 48);
    }

#elif defined(_FEC_USE_POLY_MUL2)
    // TODO: check how much performance it adds
    fec_idx_t i;
    uint32_t aa = ((uint32_t)pak_val << 0) | ((uint32_t)pak_val << 16);
    for (i = 0; i < num_y_missing; i+=2) {
        uint32_t bb = inv_arr[poly_add(pak_xy, missing_y[i])];
        if (i != num_y_missing - 1) {
            bb |= ((uint32_t)inv_arr[poly_add(pak_xy, missing_y[i+1])]) << 16;
        }
        uint32_t cc = poly_mul2_2(aa, bb);
        *((uint32_t*)&recovered[i]) ^= cc;
    }

#elif defined(_FEC_USE_POLY_MUL_CLMUL2)
    fec_idx_t i;
    for (i = 0; i < num_y_missing-1; i+=2) {
        uint32_t cc = poly_mul2(pak_val, pak_val, inv_arr[poly_add(pak_xy, missing_y[i])], inv_arr[poly_add(pak_xy, missing_y[i+1])]);
        *((uint32_t*)&recovered[i]) ^= cc;
    }
    if (i != num_y_missing) {
        recovered[i] ^= poly_mul(pak_val, inv_arr[poly_add(pak_xy, missing_y[i])]);
    }
#else
#error all cases should be covered
#endif
*/
    
}


#ifdef FEC_HAS_CLMUL64
static void PERF_DEBUG_ATTRS __fec_rx_col_op_4(uint32_t* restrict recovered, const fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_int_t* restrict inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {
    fec_idx_t i;

    _poly_t _b = _POLY_INT_TO_POLY(pak_val);
    
    for (i = 0; i < num_y_missing - 1; i+=2) {
        // uint32_t h = ((uint32_t*)missing_y)[i/2];
        // uint16_t a1 = inv_arr[poly_add(pak_xy, h)];
        // uint16_t a2 = inv_arr[poly_add(pak_xy, h >> 16)];
        uint16_t a1 = inv_arr[poly_add(pak_xy, missing_y[i])];
        uint16_t a2 = inv_arr[poly_add(pak_xy, missing_y[i+1])];

        _poly_t _a = _POLY_TWO_INT_TO_POLY(a1, a2);

        _poly_t _c = _POLY_SINGLE_MUL(_a, _b);
        
        recovered[i] ^= _POLY_GET_LOWER_32(_c);
        recovered[i+1] ^= _POLY_GET_UPPER_32(_c);
    }

    if (i == num_y_missing - 1) {
        uint16_t a1 = inv_arr[poly_add(pak_xy, missing_y[i])];
        _poly_t _a = _POLY_INT_TO_POLY(a1);
        _poly_t _c = _POLY_SINGLE_MUL(_a, _b);
        recovered[i] ^= _POLY_GET_LOWER_32(_c);
    }
}
#endif

static void PERF_DEBUG_ATTRS __fec_rx_col_op_5(u16x16* restrict recovered, const fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_int_t* restrict inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {
    fec_idx_t i;

    u16x16 shifts = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    u16x16 _b = (u16x16)(((i16x16)(pak_val << shifts)) >> 15);

    for (i = 0; i < num_y_missing; i++) {
        uint16_t a = inv_arr[poly_add(pak_xy, missing_y[i])];
        u16x16 _c = a & _b;
        recovered[i] ^= _c;
    }
}

#ifdef FEC_HAS_64_INT_VEC

#if (defined(__x86_64__) || defined(__i386__)) && defined(__MMX__)
void PERF_DEBUG_ATTRS fec_rx_col_op_mmx(__m64* restrict recovered, const fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_int_t* restrict inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {
    fec_idx_t i;

    int j;
    __m64 _pak_val[4];
    const __m64 __pak_val = _mm_set1_pi16(pak_val);
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        _pak_val[j] = _mm_srai_pi16(_mm_mullo_pi16(__pak_val, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}), 15);
    }

    for (i = 0; i < num_y_missing; i++) {
        __m64 a = _mm_set1_pi16(inv_arr[poly_add(pak_xy, missing_y[i])]);
        #pragma GCC unroll 4
        for (j = 0; j < 4; j++) {

            // recovered[i*4 + j] = _mm_xor_si64(recovered[i*4 + j], _mm_and_si64(_pak_val[j], a));

            __m64 res = _mm_and_si64(_pak_val[j], a);

            asm inline(
                "pxor %0, %1;"
                "movq %1, %0;"
            : "+m"(recovered[i*4 + j]), "+y"(res)
            : 
            : );
        }
    }
}
#endif

#if (defined(__x86_64__) || defined(__i386__)) && defined(__SSE__)
void PERF_DEBUG_ATTRS __attribute__((visibility("default"))) fec_rx_col_op_sse(__m128* restrict recovered, const fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_int_t* restrict inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {
    fec_idx_t i;

    int j;
    union {
        __m64 mmx[4];
        __m128 sse[2];
    } _pak_val_tmp;
    const __m64 __pak_val = _mm_set1_pi16(pak_val);
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        __m64 res = _mm_srai_pi16(_mm_mullo_pi16(__pak_val, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}), 15);
        asm inline(
                "movq %1, %0;"
            : "=m"(_pak_val_tmp.mmx[j])
            : "y"(res)
            : );
    }

     __m128 _pak_val_low;
     __m128 _pak_val_high;

    asm inline(
                "movaps %1, %0;"
            : "=x"(_pak_val_low)
            : "m"(_pak_val_tmp.sse[0])
            : );
    
    asm inline(
                "movaps %1, %0;"
            : "=x"(_pak_val_high)
            : "m"(_pak_val_tmp.sse[1])
            : );

    for (i = 0; i < num_y_missing; i++) {

        uint16_t a = inv_arr[poly_add(pak_xy, missing_y[i])];

        uint32_t a_2 =  (a | (((uint32_t)a) << 16));
        __m128 a_8;

        __m128 res;

        asm inline(
                "movss %[a_2], %[a_8];"
                "shufps $0, %[a_8], %[a_8];"
            : [a_8] "=x"(a_8)
            : [a_2] "m"(a_2)
            : );

        res = a_8;

        asm inline(
                "andps %[_pak_val_low], %[res];"
            : [res] "+x"(res)
            : [_pak_val_low] "x"(_pak_val_low)
            : );

        asm inline(
                "xorps %0, %1;"
                "movaps %1, %0;"
            : "+m"(recovered[2*i]), "+x"(res)
            : 
            : );

        res = a_8;

        asm inline(
                "andps %[_pak_val_high], %[res];"
            : [res] "+x"(res)
            : [_pak_val_high] "x"(_pak_val_high)
            : );

        asm inline(
                "xorps %0, %1;"
                "movaps %1, %0;"
            : "+m"(recovered[2*i + 1]), "+x"(res)
            : 
            : );

        // recovered[2*i] = _mm_xor_ps(recovered[2*i], _mm_and_ps(_pak_val_low, a_8));
        // recovered[2*i + 1] = _mm_xor_ps(recovered[2*i + 1], _mm_and_ps(_pak_val_high, a_8));
    }
}
#endif

#endif

#if defined(FEC_HAS_64BIT)
void PERF_DEBUG_ATTRS fec_rx_col_op_reg64(uint64_t* restrict recovered, const fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_int_t* restrict inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {

    int j;
    uint64_t _pak_val[4];

    fec_idx_t i;
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        _pak_val[j] = pak_val * ((1ULL<<((16+1)*(4*j+0) - 64*j)) + (1ULL<<((16+1)*(4*j+1) - 64*j)) + (1ULL<<((16+1)*(4*j+2) - 64*j)) + (1ULL<<((16+1)*(4*j+3) - 64*j)));
        _pak_val[j] >>= (16-1);
        _pak_val[j] &= (1ULL<<(16*0)) | (1ULL<<(16*1)) | (1ULL<<(16*2)) | (1ULL<<(16*3));
        _pak_val[j] *= (1<<16) - 1;
    }

    for (i = 0; i < num_y_missing; i++) {
        uint64_t b = (inv_arr[poly_add(pak_xy, missing_y[i])] * ((1ULL<<(16*0)) + (1ULL<<(16*1)) + (1ULL<<(16*2)) + (1ULL<<(16*3))));

        #pragma GCC unroll 4
        for (j = 0; j < 4; j++) {
            recovered[i*4+j] ^= _pak_val[j] & b;
        }
    }
}
#endif

#if defined(FEC_HAS_32BIT)
void PERF_DEBUG_ATTRS fec_rx_col_op_reg32(uint32_t* restrict recovered, const fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_int_t* restrict inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {

    int j;
    uint32_t _pak_val[8];

    fec_idx_t i;
    
    #pragma GCC unroll 8
    for (j = 0; j < 8; j++) {
        _pak_val[j] = ((uint16_t)(((int16_t)(pak_val << (16 - 1 - (j*2 + 0)))) >> 15)) | (((uint32_t)((uint16_t)(((int16_t)(pak_val << (16 - 1 - (j*2 + 1)))) >> 15))) << 16);
    }

    for (i = 0; i < num_y_missing; i++) {
        uint16_t b = inv_arr[poly_add(pak_xy, missing_y[i])];
        uint32_t b_2 = b | (((uint32_t)b) << 16);

        #pragma GCC unroll 8
        for (j = 0; j < 8; j++) {
            recovered[i*8+j] ^= _pak_val[j] & b_2;
        }
    }
}
#endif

fec_int_t PERF_DEBUG_ATTRS __fec_rx_row_op(const unaligend_fec_int_t * const* pak_arr, fec_idx_t num, size_t ii, const fec_int_t* xy_arr, fec_int_t missing_y_i, const fec_int_t *inv_arr, fec_int_t ones_pak_ii) {
    fec_idx_t j;
    //fec_int_t res = ones_pak_ii;

    inv_arr -= 1;

    // for (j = 0; j < num; j++) {
    //     res = poly_add(res, poly_mul(pak_arr[j][ii], inv_arr[poly_add(xy_arr[j], missing_y_i)]));
    // }

    u32x8 sum1 = {0};

    u32x8 cc1 = {0,1,2,3,4,5,6,7};
    u32x8 cc2 = {8,9,10,11,12,13,14,15};

    for (j = 0; j < num; j++) {
        fec_int_t a = pak_arr[j][ii];
        fec_int_t b = inv_arr[poly_add(xy_arr[j], missing_y_i)];

        u32x8 aa1 = {a,a,a,a,a,a,a,a};
        u32x8 aa2 = {a,a,a,a,a,a,a,a};
        u32x8 bb1 = {b,b,b,b,b,b,b,b};
        u32x8 bb2 = {b,b,b,b,b,b,b,b};
        
        aa1 <<= cc1;
        bb1 = -((bb1 >> cc1) & 1);
        aa2 <<= cc2;
        bb2 = -((bb2 >> cc2) & 1);
        
        sum1 ^= (aa1 & bb1) ^ (aa2 & bb2);
    }

    u32x4 sum2 = (u32x4)(my_mm256_extracti128_si256((__m256i)sum1, 0) ^ my_mm256_extracti128_si256((__m256i)sum1, 1));
    uint32_t sum = sum2[0] ^ sum2[1] ^ sum2[2] ^ sum2[3];

    uint16_t res = sum;
    uint32_t carry = sum >> 16;

    carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

    res ^= carry;

    carry = carry >> 16;
    
    res ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

    return res ^ ones_pak_ii;

    // //fec_idx_t i;
    // u16x16 res_arr = {0};
    // for (j = 0; j < num; j+=16) {
    //     u16x16 aa,bb;

    //     // fec_idx_t j;
    //     // fec_idx_t num_to_copy = MIN(16, num_y_missing - i);
    //     // for(j = 0; j < num_to_copy; j++) {
    //     //     bb[j] = inv_arr[poly_add(pak_xy, missing_y[i+j])];
    //     // }

    //     uint16_t __attribute__((aligned(32))) tmp_arr[16] = {0};
    //     uint16_t __attribute__((aligned(32))) tmp_arr2[16] = {0};
    //     fec_idx_t jj;
    //     fec_idx_t num_to_copy = MIN(16U, num - j);
    //     for(jj = 0; jj < num_to_copy; jj++) {
    //         tmp_arr[jj] = inv_arr[poly_add(xy_arr[j+jj], missing_y_i)];
    //         tmp_arr2[jj] = pak_arr[j+jj][ii];
    //     }
    //     aa = (u16x16)_mm256_load_si256((__m256i*)tmp_arr);
    //     bb = (u16x16)_mm256_load_si256((__m256i*)tmp_arr2);
        

    //     // bb[0] = inv_arr[poly_add(pak_xy, missing_y[i+0])];
    //     // bb[1] = inv_arr[poly_add(pak_xy, missing_y[i+1])];
    //     // bb[2] = inv_arr[poly_add(pak_xy, missing_y[i+2])];
    //     // bb[3] = inv_arr[poly_add(pak_xy, missing_y[i+3])];
    //     // bb[4] = inv_arr[poly_add(pak_xy, missing_y[i+4])];
    //     // bb[5] = inv_arr[poly_add(pak_xy, missing_y[i+5])];
    //     // bb[6] = inv_arr[poly_add(pak_xy, missing_y[i+6])];
    //     // bb[7] = inv_arr[poly_add(pak_xy, missing_y[i+7])];
    //     // bb[8] = inv_arr[poly_add(pak_xy, missing_y[i+8])];
    //     // bb[9] = inv_arr[poly_add(pak_xy, missing_y[i+9])];
    //     // bb[10] = inv_arr[poly_add(pak_xy, missing_y[i+10])];
    //     // bb[11] = inv_arr[poly_add(pak_xy, missing_y[i+11])];
    //     // bb[12] = inv_arr[poly_add(pak_xy, missing_y[i+12])];
    //     // bb[13] = inv_arr[poly_add(pak_xy, missing_y[i+13])];
    //     // bb[14] = inv_arr[poly_add(pak_xy, missing_y[i+14])];
    //     // bb[15] = inv_arr[poly_add(pak_xy, missing_y[i+15])];
    //     u16x16 cc = poly_mul16(aa, bb);

    //     res_arr ^= cc;

    //     // recovered[i+0] ^= cc[0];
    //     // recovered[i+1] ^= cc[1];
    //     // recovered[i+2] ^= cc[2];
    //     // recovered[i+3] ^= cc[3];
    //     // recovered[i+4] ^= cc[4];
    //     // recovered[i+5] ^= cc[5];
    //     // recovered[i+6] ^= cc[6];
    //     // recovered[i+7] ^= cc[7];
    //     // recovered[i+8] ^= cc[8];
    //     // recovered[i+9] ^= cc[9];
    //     // recovered[i+10] ^= cc[10];
    //     // recovered[i+11] ^= cc[11];
    //     // recovered[i+12] ^= cc[12];
    //     // recovered[i+13] ^= cc[13];
    //     // recovered[i+14] ^= cc[14];
    //     // recovered[i+15] ^= cc[15];
    // }

    // fec_idx_t res = ones_pak_ii;
    // for (j = 0; j < 16; j++) {
    //     res ^= res_arr[j];
    // }

    // return res;



    // __m128i tmp1 = _mm256_extracti128_si256((__m256i)res_arr, 0);
    // __m128i tmp2 = _mm256_extracti128_si256((__m256i)res_arr, 1);
    // tmp1 ^= tmp2;

    // uint64_t tmp3 = _mm_extract_epi64(tmp1, 0);
    // uint64_t tmp4 = _mm_extract_epi64(tmp1, 1);
    // tmp3 ^= tmp4;

    // tmp3 ^= (tmp3>>32);

    // tmp3 ^= (tmp3>>16);

    // return tmp3 ^ ones_pak_ii;
}



#ifdef PERF_DEBUG
#ifdef _WIN32
#include <realtimeapiset.h>
#include <processthreadsapi.h>
static uint64_t get_timestamp() {
    // struct timespec tp = {0};
    // //CHECK(clock_gettime(CLOCK_MONOTONIC, &tp) == 0);
    // clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp);
    // return (tp.tv_sec*1000000000ULL) + tp.tv_nsec;
    uint64_t t;
    QueryThreadCycleTime(GetCurrentThread(), &t);
    t /= 4;
    return t;
}
#else
static uint64_t get_timestamp() {
    struct timespec tp = {0};
    //CHECK(clock_gettime(CLOCK_MONOTONIC, &tp) == 0);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp);
    return (tp.tv_sec*1000000000ULL) + tp.tv_nsec;
}
#endif
#endif

#ifdef FEC_USER_GIVEN_BUFFER
static void PERF_DEBUG_ATTRS memswap(unaligend_fec_int_t* a, unaligend_fec_int_t* b, size_t size) {
    // TODO: i am hoping compiler will optimize this

    size_t i;
    for (i = 0; i < size; i++) {
        fec_int_t tmp = a[i];
        a[i] = b[i];
        b[i] = tmp;
    }
}
#endif

char gggg[] = "aaaaaaaaa";
char hhhh[] = "bbbbbbbbb";

#include "preprocessor.h"

void PERF_DEBUG_ATTRS __fec_rx_one_op3(fec_idx_t num, fec_int_t* restrict present_y, unaligend_fec_int_t * restrict * restrict pak_arr, fec_perf_int_t* tmp_recovered_ints, fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_inv_cache_t* restrict inv_cache, size_t ii) {
    fec_idx_t j;
    for (j = 0; j < num; j++) {
        fec_int_t xy_j = present_y[j];
        fec_int_t xy_pak_arr_j_ii = pak_arr[j][ii];
        // for (i = 0; i < num_y_missing; i++) {
        //     fec_int_t missing_y_i = missing_y[i];
        //     tmp_recovered_ints[i] ^= poly_mul(xy_pak_arr_j_ii, _fec_inv(inv_cache, poly_add(xy_j, missing_y_i)));
        // }
        __fec_rx_col_op(tmp_recovered_ints, missing_y, num_y_missing, inv_cache->inv_arr, xy_pak_arr_j_ii, xy_j);
    }

    /*fec_idx_t i;
    for (i = 0; i < num_y_missing; i++) {
        uint16_t res = tmp_recovered_ints[i];
        uint32_t carry = tmp_recovered_ints[i] >> 16;

        carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

        res ^= carry;

        carry = carry >> 16;
        
        res ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);
        
        tmp_recovered_ints[i] = res;
    }*/
    // fec_idx_t i;
    // u32x8 shifts1 = {0, 2, 4, 6, 8, 10, 12, 14};
    // u32x8 shifts2 = {1, 3, 5, 7, 9, 11, 13, 15};
    // u64x2 _poly = {POLY_G, 0};

    // for (i = 0; i < num_y_missing; i++) {
    //     u16x16 res = tmp_recovered_ints[i];
    //     u32x8 res1 = (((u32x8)res) & 0xFFFF) << shifts1;
    //     u32x8 res2 = (((u32x8)res) >> 16) << shifts2;

    //     // u64x2 res1_1 = (u64x2)(_mm256_extracti128_si256((__m256i)res1, 0) ^ _mm256_extracti128_si256((__m256i)res1, 1));
    //     // u64x2 res2_1 = (u64x2)(_mm256_extracti128_si256((__m256i)res2, 0) ^ _mm256_extracti128_si256((__m256i)res2, 1));

    //     // uint64_t res1_2 = res1_1[0] ^ res1_1[1];
    //     // uint64_t res2_2 = res2_1[0] ^ res2_1[1];

    //     // uint32_t res1_x = (res1_2 >> 32) ^ res1_2;
    //     // uint32_t res2_x = (res2_2 >> 32) ^ res2_2;

    //     uint32_t res1_x = res1[0] ^ res1[1] ^ res1[2] ^ res1[3] ^ res1[4] ^ res1[5] ^ res1[6] ^ res1[7];
    //     uint32_t res2_x = res2[0] ^ res2[1] ^ res2[2] ^ res2[3] ^ res2[4] ^ res2[5] ^ res2[6] ^ res2[7];

    //     uint32_t res3 = res1_x ^ res2_x;

        
    //     u32x4 _c = {res3};
    //     u32x4 _d = _c >> 16;
    //     _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    //     _c ^= _d;
    //     _d >>= 16;
    //     _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    //     _c ^= _d;

    //     ((uint16_t*)tmp_recovered_ints)[i] = ((u16x8)_c)[0];

    //     // uint16_t res4 = res3;
    //     // uint32_t carry = res3 >> 16;

    //     // carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

    //     // res4 ^= carry;

    //     // carry = carry >> 16;
        
    //     // res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);
        
    //     // tmp_recovered_ints[i][0] = res4;
    // }
}

void __fec_rx_one_op(fec_idx_t num, fec_int_t* present_y, unaligend_fec_int_t **pak_arr, fec_int_t *tmp_recovered_ints, fec_int_t *missing_y, fec_idx_t num_y_missing, const fec_inv_cache_t* inv_cache, size_t ii) {
    // fec_idx_t j;
    // for (j = 0; j < num; j++) {
    //     fec_int_t xy_j = present_y[j];
    //     fec_int_t xy_pak_arr_j_ii = pak_arr[j][ii];
    //     // for (i = 0; i < num_y_missing; i++) {
    //     //     fec_int_t missing_y_i = missing_y[i];
    //     //     tmp_recovered_ints[i] ^= poly_mul(xy_pak_arr_j_ii, _fec_inv(inv_cache, poly_add(xy_j, missing_y_i)));
    //     // }
    //     volatile char* volatile aaa = gggg;
    //     aaa = aaa + 1;
    //     *aaa = *aaa;
    //     __fec_rx_col_op(tmp_recovered_ints, missing_y, num_y_missing, inv_cache->inv_arr, xy_pak_arr_j_ii, xy_j);

    //     aaa = hhhh;
    //     aaa = aaa + 1;
    //     *aaa = *aaa;
    // }

    fec_idx_t i,j,k,m;

    fec_int_t prefetch[32];

    const fec_int_t* inv_arr = inv_cache->inv_arr - 1;

    u16x16 cc1 = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    u16x16 cc2 = 16 - cc1;

    for (j = 0; j < num; j += sizeof(prefetch)/sizeof(prefetch[0])) {
        fec_idx_t num_k = MIN(sizeof(prefetch)/sizeof(prefetch[0]), num - j);

        for(k = 0; k < num_k; k++) {
            prefetch[k] = pak_arr[j + k][ii];
        }

        #define NEW_MUL_NUM 8

        for(i = 0; i < num_y_missing; i += NEW_MUL_NUM) {
            
            //u16x16 sum1[NEW_MUL_NUM] = {0};
            #define MMM(i, _) u16x16 _GLUE(sum1_,i) = {0};
            EVAL(REPEAT(NEW_MUL_NUM, MMM, ~));
            #undef MMM

            fec_idx_t num_m = MIN(NEW_MUL_NUM, num_y_missing - i);

            for (k = 0; k < num_k; k++) {
                fec_int_t b = prefetch[k];
                u16x16 bb1 = {b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b};
                bb1 = -((bb1 >> cc1) & 1);
                fec_int_t pres_y = present_y[j+k];

                // #pragma GCC unroll NEW_MUL_NUM
                // for(m = 0; m < sizeof(sum1)/sizeof(sum1[0]); m++) {
                //     fec_int_t missing_y_i = missing_y[m < num_m ? i+m: i];
                //     fec_int_t a = inv_arr[poly_add(pres_y, missing_y_i)];

                //     u16x16 aa1 = {a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a};
                    
                //     sum1[m] ^= (aa1 & bb1);
                // }

                #define MMM(m2, _) { \
                    fec_int_t missing_y_i = missing_y[m2 < num_m ? i+m2: i]; \
                    fec_int_t a = inv_arr[poly_add(pres_y, missing_y_i)]; \
                    \
                    u16x16 aa1 = {a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a}; \
                    \
                    _GLUE(sum1_,m2) ^= (aa1 & bb1); \
                }
                EVAL(REPEAT(NEW_MUL_NUM, MMM, ~));
                #undef MMM
            }

            // #pragma GCC unroll NEW_MUL_NUM
            // for(m = 0; m < NEW_MUL_NUM; m++) {
            //     u16x16 sum2, sum3;
            //     sum2 = sum1[m] << cc1;
            //     sum3 = sum1[m] >> cc2;
            //     sum3[0] = 0;

            //     u16x8 sum2_2, sum3_2;
            //     sum2_2 = (u16x8)(_mm256_extracti128_si256((__m256i)sum2, 0) ^ _mm256_extracti128_si256((__m256i)sum2, 1));
            //     sum3_2 = (u16x8)(_mm256_extracti128_si256((__m256i)sum3, 0) ^ _mm256_extracti128_si256((__m256i)sum3, 1));

            //     uint64_t sum2_3, sum3_3;
            //     sum2_3 = ((u64x2)sum2_2)[0] ^ ((u64x2)sum2_2)[1];
            //     sum3_3 = ((u64x2)sum3_2)[0] ^ ((u64x2)sum3_2)[1];

            //     uint32_t sum2_4, sum3_4;
            //     sum2_4 = sum2_3 ^ (sum2_3 >> 32);
            //     sum3_4 = sum3_3 ^ (sum3_3 >> 32);

            //     uint16_t res = sum2_4 ^ (sum2_4 >> 16);
            //     uint32_t carry = (uint16_t)(sum3_4 ^ (sum3_4 >> 16));

            //     carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

            //     res ^= carry;

            //     carry = carry >> 16;
                
            //     res ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

            //     tmp_recovered_ints[i+m] ^= res;
            // }

            #define MMM(m2, _) { \
                u16x16 sum2, sum3; \
                sum2 = _GLUE(sum1_,m2) << cc1; \
                sum3 = _GLUE(sum1_,m2) >> cc2; \
                sum3[0] = 0; \
                \
                u16x8 sum2_2, sum3_2; \
                sum2_2 = (u16x8)(my_mm256_extracti128_si256((__m256i)sum2, 0) ^ my_mm256_extracti128_si256((__m256i)sum2, 1)); \
                sum3_2 = (u16x8)(my_mm256_extracti128_si256((__m256i)sum3, 0) ^ my_mm256_extracti128_si256((__m256i)sum3, 1)); \
                \
                uint64_t sum2_3, sum3_3; \
                sum2_3 = ((u64x2)sum2_2)[0] ^ ((u64x2)sum2_2)[1]; \
                sum3_3 = ((u64x2)sum3_2)[0] ^ ((u64x2)sum3_2)[1]; \
                \
                uint32_t sum2_4, sum3_4; \
                sum2_4 = sum2_3 ^ (sum2_3 >> 32); \
                sum3_4 = sum3_3 ^ (sum3_3 >> 32); \
                \
                uint16_t res = sum2_4 ^ (sum2_4 >> 16); \
                uint32_t carry = (uint16_t)(sum3_4 ^ (sum3_4 >> 16)); \
                \
                carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5); \
                \
                res ^= carry; \
                \
                carry = carry >> 16; \
                \
                res ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5); \
                \
                tmp_recovered_ints[i+m2] ^= res; \
            }
            EVAL(REPEAT(NEW_MUL_NUM, MMM, ~));
            #undef MMM
        }
    }

}

#define UNROLL_1(x) x
#define UNROLL_2(x) UNROLL_1(x) UNROLL_1(x)
#define UNROLL_4(x) UNROLL_2(x) UNROLL_2(x)
#define UNROLL_8(x) UNROLL_4(x) UNROLL_4(x)
#define UNROLL_16(x) UNROLL_8(x) UNROLL_8(x)

#define UNROLL_N(N, x) \
    if(N & 16) { UNROLL_16(x) } \
    if(N & 8) { UNROLL_8(x) } \
    if(N & 4) { UNROLL_4(x) } \
    if(N & 2) { UNROLL_2(x) } \
    if(N & 1) { UNROLL_1(x) }


#define ONEREP_1 a
#define ONEREP_2 ONEREP_1##ONEREP_1
#define ONEREP_4 ONEREP_2##ONEREP_2
#define ONEREP_8 ONEREP_4##ONEREP_4
#define ONEREP_16 ONEREP_8##ONEREP_8


#define UNROLL_LOOP(N, init, inc, x) init UNROLL_N(N, {x inc})

#ifdef __PCLMUL__

void __attribute__((noinline)) __fec_rx_one_op2(fec_idx_t num, fec_int_t* present_y, unaligend_fec_int_t **pak_arr, fec_int_t *tmp_recovered_ints, fec_int_t *missing_y, fec_idx_t num_y_missing, const fec_inv_cache_t* inv_cache, size_t ii) {

    fec_idx_t i,j,k,m;

    fec_int_t prefetch[32];

    u64x2 poly_g = {POLY_G, 0};

    const fec_int_t* inv_arr = inv_cache->inv_arr - 1;

    #pragma GCC unroll 1
    for (j = 0; j < num; j += sizeof(prefetch)/sizeof(prefetch[0])) {
        fec_idx_t num_k = MIN(sizeof(prefetch)/sizeof(prefetch[0]), num - j);

        #pragma GCC unroll 1
        for(k = 0; k < num_k; k++) {
            prefetch[k] = pak_arr[j + k][ii];
        }

        #define NUM_UNROLLS 8

        #pragma GCC unroll 1
        for(i = 0; i < num_y_missing; i += 2*NUM_UNROLLS) {
            
            uint64_t sum1[NUM_UNROLLS] = {0};
            #define MMM(i, _) uint64_t _GLUE(sum1_,i) = 0;
            EVAL(REPEAT(NUM_UNROLLS, MMM, ~));
            #undef MMM

            fec_idx_t num_m = MIN(2*NUM_UNROLLS, num_y_missing - i);

            fec_int_t missing_y_prefetch[2*NUM_UNROLLS];

            //memcpy(missing_y_prefetch, &missing_y[i], num_m * sizeof(fec_int_t));
            for (m = 0; m < num_m; m++) {
                missing_y_prefetch[m] = missing_y[i+m];
            }
            for (m = num_m; m < 2*NUM_UNROLLS; m++) {
                missing_y_prefetch[m] = missing_y[i];
            }
            //memset(&missing_y_prefetch[num_m], 0, (2*NUM_UNROLLS - num_m)*sizeof(fec_int_t));

            #pragma GCC unroll 1
            for (k = 0; k < num_k; k++) {
                fec_int_t b = prefetch[k];
                u64x2 bb = {0};
                bb[0] = b;
                fec_int_t pres_y = present_y[j+k];

                // #pragma unroll NUM_UNROLLS
                // for(m = 0; m < 2*sizeof(sum1)/sizeof(sum1[0]); m+=2) {
                //     fec_int_t missing_y_i = missing_y[m < num_m ? i+m:i];
                //     fec_int_t missing_y_i_2 = missing_y[m + 1 < num_m ? i+m+1:i];
                //     fec_int_t a = inv_arr[poly_add(present_y[j+k], missing_y_i)];
                //     fec_int_t a2 = inv_arr[poly_add(present_y[j+k], missing_y_i_2)];
                //     u32x4 aa = {0};
                    
                //     aa[0] = a;
                //     aa[1] = a2;

                //     sum1[m/2] ^= ((u64x2)_mm_clmulepi64_si128((__m128i)aa, (__m128i)bb, 0))[0];
                // }

                // UNROLL_LOOP(NUM_UNROLLS, {m = 0;}, {m += 2;}, {
                //     fec_int_t missing_y_i = missing_y[m < num_m ? i+m:i];
                //     fec_int_t missing_y_i_2 = missing_y[m + 1 < num_m ? i+m+1:i];
                //     fec_int_t a = inv_arr[poly_add(present_y[j+k], missing_y_i)];
                //     fec_int_t a2 = inv_arr[poly_add(present_y[j+k], missing_y_i_2)];
                //     u32x4 aa = {0};
                    
                //     aa[0] = a;
                //     aa[1] = a2;

                //     sum1[m/2] ^= ((u64x2)_mm_clmulepi64_si128((__m128i)aa, (__m128i)bb, 0))[0];
                // });

                #define MMM(m2, _) { \
                    fec_int_t missing_y_i = missing_y_prefetch[m2*2 < num_m ? m2*2:0]; \
                    fec_int_t missing_y_i_2 = missing_y_prefetch[m2*2+1 < num_m ? m2*2+1:0]; \
                    fec_int_t a = inv_arr[poly_add(pres_y, missing_y_i)]; \
                    fec_int_t a2 = inv_arr[poly_add(pres_y, missing_y_i_2)]; \
                    u32x4 aa = {0}; \
                     \
                    aa[0] = a; \
                    aa[1] = a2; \
                     \
                    _GLUE(sum1_,m2) ^= ((u64x2)_mm_clmulepi64_si128((__m128i)aa, (__m128i)bb, 0))[0]; \
                }
                EVAL(REPEAT(NUM_UNROLLS, MMM, ~));
                #undef MMM
                

            }

            
            //poly_g[0] = POLY_G;

            // #pragma unroll NUM_UNROLLS
            // for(m = 0; m < num_m; m+=2) {
            //     uint64_t sum = sum1[m/2];
            //     u64x2 aa = {0};
            //     aa[0] = sum;

            //     u32x4 bb = (u32x4)_mm_clmulepi64_si128((__m128i)(((u32x4)aa) >> 16), (__m128i)poly_g, 0);
            //     bb ^= (u32x4)_mm_clmulepi64_si128((__m128i)(bb >> 16), (__m128i)poly_g, 0);
            //     bb ^= (u32x4)aa;

            //     uint16_t res1 = ((u16x8)bb)[0];
            //     uint16_t res2 = ((u16x8)bb)[2];

            //     tmp_recovered_ints[i+m] ^= res1;
            //     tmp_recovered_ints[i+m+1] ^= res2;
            // }

            // UNROLL_LOOP(NUM_UNROLLS, {m = 0;}, {m += 2;}, {
            //     if (m < num_m) {
            //         uint64_t sum = sum1[m/2];
            //         u64x2 aa = {0};
            //         aa[0] = sum;

            //         u32x4 bb = (u32x4)_mm_clmulepi64_si128((__m128i)(((u32x4)aa) >> 16), (__m128i)poly_g, 0);
            //         bb ^= (u32x4)_mm_clmulepi64_si128((__m128i)(bb >> 16), (__m128i)poly_g, 0);
            //         bb ^= (u32x4)aa;

            //         uint16_t res1 = ((u16x8)bb)[0];
            //         uint16_t res2 = ((u16x8)bb)[2];

            //         tmp_recovered_ints[i+m] ^= res1;
            //         tmp_recovered_ints[i+m+1] ^= res2;
            //     }
            // });

            #define MMM(m2, _) { \
                /*if (m2*2 < num_m) {*/ \
                    uint64_t sum = _GLUE(sum1_,m2); \
                    u64x2 aa = {0}; \
                    aa[0] = sum; \
                     \
                    u32x4 bb = (u32x4)_mm_clmulepi64_si128((__m128i)(((u32x4)aa) >> 16), (__m128i)poly_g, 0); \
                    bb ^= (u32x4)_mm_clmulepi64_si128((__m128i)(bb >> 16), (__m128i)poly_g, 0); \
                    bb ^= (u32x4)aa; \
                     \
                    uint16_t res1 = ((u16x8)bb)[0]; \
                    uint16_t res2 = ((u16x8)bb)[2]; \
                     \
                    tmp_recovered_ints[i+m2*2] ^= res1; \
                    tmp_recovered_ints[i+m2*2+1] ^= res2; \
                /*}*/ \
            }
            EVAL(REPEAT(NUM_UNROLLS, MMM, ~));
            #undef MMM
            // #define MMM(m2, _) { \
            //     _GLUE(the_label, m2): \
            //     case m2*2+1: \
            //     case m2*2+2: { \
            //         uint64_t sum = _GLUE(sum1_,m2); \
            //         u64x2 aa = {0}; \
            //         aa[0] = sum; \
            //          \
            //         u32x4 bb = (u32x4)_mm_clmulepi64_si128((__m128i)(((u32x4)aa) >> 16), (__m128i)poly_g, 0); \
            //         bb ^= (u32x4)_mm_clmulepi64_si128((__m128i)(bb >> 16), (__m128i)poly_g, 0); \
            //         bb ^= (u32x4)aa; \
            //          \
            //         uint16_t res1 = ((u16x8)bb)[0]; \
            //         uint16_t res2 = ((u16x8)bb)[2]; \
            //          \
            //         tmp_recovered_ints[i+m2*2] ^= res1; \
            //         tmp_recovered_ints[i+m2*2+1] ^= res2; \
            //         IF(m2)(goto GLUE(the_label, DEC(m2));, break;) \
            //     } \
            // }
            // switch(num_m) {
            // EVAL(REPEAT(NUM_UNROLLS, MMM, ~));
            // }
            // #undef MMM
        }
    }

}

void __attribute__((noinline)) __attribute__((visibility("default"))) blabla(fec_idx_t num_y_missing, fec_idx_t n, size_t pak_len, bool has_one_row) {
    size_t i;
    uint64_t gggg[16] = {0};
    for(i = 0; i < num_y_missing*(n - has_one_row)*pak_len; i++) {
        u64x2 aa = {0};
        u64x2 bb = {0};
        u64x2 cc;
        aa[0] = i << 1;
        bb[0] = i << 2;
        cc = (u64x2)_mm_clmulepi64_si128((__m128i)aa, (__m128i)bb, 0);
        gggg[i & 0xF] ^= cc[0];
    }
    volatile uint64_t aaaa = 0;
    for(i = 0; i < 16; i++) {
        aaaa += gggg[i];
    }
    aaaa = aaaa;
}

#endif

#ifndef FEC_USER_GIVEN_BUFFER
#define _PAK_ARR_PARAMS unaligend_fec_int_t* restrict * restrict y_pak_arr, unaligend_fec_int_t* restrict * restrict x_pak_arr
#else
#define _PAK_ARR_PARAMS unaligend_fec_int_t* restrict y_paks_buf, unaligend_fec_int_t* restrict x_pak_arr size_t pak_len
#endif

void PERF_DEBUG_ATTRS __fec_rx_one_op4(fec_idx_t num_xy, fec_int_t* restrict present_y, _PAK_ARR_PARAMS, uint32_t* restrict tmp_recovered_ints, fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_inv_cache_t* restrict inv_cache, size_t out_pak_idx, fec_int_t ones_pak_val);
void PERF_DEBUG_ATTRS __fec_rx_one_op5(fec_idx_t num_xy, fec_int_t* restrict present_y, _PAK_ARR_PARAMS, u16x16* restrict tmp_recovered_ints, fec_int_t* restrict missing_y, fec_idx_t num_y_missing, const fec_inv_cache_t* restrict inv_cache, size_t out_pak_idx, fec_int_t ones_pak_val);

fec_status_t fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state, const fec_inv_cache_t *inv_cache) {
    fec_idx_t n = rx_state->n;
    size_t pak_len = rx_state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present;
    fec_int_t *present_x;
    fec_int_t *present_y;
    fec_int_t *missing_y;
    unaligend_fec_int_t *ones_pak = NULL;
    fec_idx_t num_y_present = rx_state->num_info;

    fec_idx_t i, j;
    fec_idx_t y_i;
    size_t ii;

    if (!_fec_can_recover(rx_state)) {
        return FEC_STATUS_MORE_PACKETS_NEEDED;
    }

    if (num_y_missing == 0) {
        goto reorder_packets;
    }

    // n == 0 -> num_y_missing == 0
    // k == 0 -> no redundants and can recover -> num_y_missing == 0
    // here n,k > 0
    // max_x = n + k - 2 -> n + k - 1 = max_x + 1
    // inv_cache->inv_arr_sz < n + k - 1 <-> inv_cache->inv_arr_sz <= max_x
    if (inv_cache->inv_arr_sz <= rx_state->max_x) {
        return FEC_STATUS_INV_CACHE_TOO_SMALL;
    }

#ifndef FEC_USER_GIVEN_BUFFER
    has_one_row = (rx_state->ones_pak != NULL);
    num_x_present = num_y_missing - has_one_row;
    if (has_one_row) {
        if(num_x_present > 0) {
            rx_state->pak_arr[num_y_present] = rx_state->pak_arr[n - 1];
            rx_state->pak_xy_arr[num_y_present] = rx_state->pak_xy_arr[n - 1];
        }
        rx_state->pak_arr[n - 1] = rx_state->ones_pak;
    }

#else
    has_one_row = rx_state->has_one_pak;
    num_x_present = num_y_missing - has_one_row;
    if (has_one_row && rx_state->ones_pak_idx != 0) {
        // we need to move the ones pak
        memswap(&rx_state->pak_buffer[(n - 1 - rx_state->ones_pak_idx)*pak_len], &rx_state->pak_buffer[(n - 1)*pak_len], pak_len);
        rx_state->pak_xy_arr[n - 1 - rx_state->ones_pak_idx] = rx_state->pak_xy_arr[n - 1];
        //rx_state->ones_pak_idx = 0;
    }
#endif

    present_x = &rx_state->pak_xy_arr[num_y_present];
    present_y = rx_state->pak_xy_arr;
    missing_y = rx_state->missing_y;

#ifndef FEC_USER_GIVEN_BUFFER
    unaligend_fec_int_t** x_pak_arr = &rx_state->pak_arr[num_y_present];
    unaligend_fec_int_t** y_pak_arr = rx_state->pak_arr;
    ones_pak = rx_state->ones_pak;
#define _GET_X_PAK(idx) (x_pak_arr[(idx)])
#define _GET_Y_PAK(idx) (y_pak_arr[(idx)])
#define _GET_XY_PAK(idx) _GET_Y_PAK(idx)
#define _PAK_ARR_PARAMS unaligend_fec_int_t* restrict * restrict y_pak_arr, unaligend_fec_int_t* restrict * restrict x_pak_arr
#else
    unaligend_fec_int_t* x_paks_buf = &rx_state->pak_buffer[num_y_present * pak_len];
    unaligend_fec_int_t* y_paks_buf = rx_state->pak_buffer;
    ones_pak = &rx_state->pak_buffer[(n - 1)*pak_len];
#define _GET_X_PAK(idx) (&x_paks_buf[(idx)*pak_len])
#define _GET_Y_PAK(idx) (&y_paks_buf[(idx)*pak_len])
#define _GET_XY_PAK(idx) _GET_Y_PAK(idx)
#define _PAK_ARR_PARAMS unaligend_fec_int_t* restrict y_paks_buf, unaligend_fec_int_t* restrict x_pak_arr size_t pak_len
#endif

#ifdef PERF_DEBUG
    uint64_t start_time, end_time;
    start_time = get_timestamp();
#endif

    for (y_i = 0, i = 0; i < num_y_missing; y_i++) {
        if ((rx_state->received_paks_bitmap[y_i / 8] & (1<<(y_i & (8-1)))) == 0) {
            missing_y[i] = y_i;
            i++;
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.1--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    for (i = 0; i < num_y_present; i++) {
        y_i = present_y[i];

        fec_int_t pi_ycomp_y_div_ycomp_x_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, poly_add(y_i, missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, _fec_inv(inv_cache, poly_add(y_i, present_x[j])));
        }

        unaligend_fec_int_t* pak = _GET_Y_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_ycomp_y_div_ycomp_x_i);
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.2--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    for (i = 0; i < num_x_present; i++) {
        fec_int_t pi_xy_div_xx_i = 1;
        fec_int_t x_i = present_x[i];
        for (j = 0; j < num_y_missing; j++) {
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, poly_add(x_i, missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, _fec_inv(inv_cache, poly_add(x_i, present_x[j])));
        }

        unaligend_fec_int_t *pak = _GET_X_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_xy_div_xx_i);
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.3--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    fec_perf_int_t* tmp_recovered_ints = (fec_perf_int_t*)(((uintptr_t)rx_state->tmp_recovered_ints + __alignof__(fec_perf_int_t) - 1) & (-__alignof__(fec_perf_int_t)));
    const fec_int_t* inv_arr = inv_cache->inv_arr - 1;

    for (ii = 0; ii < pak_len; ii++) {
        

        fec_int_t ones_pak_ii = 0;
        if (has_one_row) {
            ones_pak_ii = ones_pak[ii];
        }

        // for (i = 0; i < num_y_missing; i++) {
        //     fec_int_t missing_y_i = missing_y[i];
        //     fec_int_t res;

        //     // res = ones_pak_ii;
        //     // for (j = 0; j < num_y_present + num_x_present; j++) {
        //     //     res = poly_add(res, poly_mul(_GET_XY_PAK(j)[ii], _fec_inv(inv_cache, poly_add(present_y[j], missing_y_i))));
        //     // }

        //     res = __fec_rx_row_op((const unaligend_fec_int_t* const*)y_pak_arr, num_y_present + num_x_present, ii, present_y, missing_y_i, inv_cache->inv_arr, ones_pak_ii);

        //     tmp_recovered_ints[i] = res;
        // }

        // for (i = 0; i < num_y_missing; i++) {
        //     _GET_X_PAK(i)[ii] = ((fec_int_t*)tmp_recovered_ints)[i];
        // }

#if defined(FEC_HAS_CLMUL64)
        {
            for (i = 0; i < num_y_missing; i++) {
                tmp_recovered_ints[i] = ones_pak_ii;
            }

            for (j = 0; j < num_y_present + num_x_present; j++) {
                fec_int_t xy_j = present_y[j];
                fec_int_t xy_pak_arr_j_val = _GET_XY_PAK(j)[ii];
                __fec_rx_col_op_4(tmp_recovered_ints, missing_y, num_y_missing, inv_arr, xy_pak_arr_j_val, xy_j);
            }

            _poly_t _poly = _POLY_INT_TO_POLY(POLY_G);

            for (i = 0; i < num_y_missing - 1; i += 2) {

                _poly_t _c = _POLY_TWO_INT_TO_POLY(tmp_recovered_ints[i], tmp_recovered_ints[i+1]);
                _poly_t _d = _POLY_TWO_SHIFT_RIGHT(_c, 16);
                _d = _POLY_SINGLE_MUL(_d, _poly);
                _c ^= _d;
                _d = _POLY_TWO_SHIFT_RIGHT(_d, 16);
                _d = _POLY_SINGLE_MUL(_d, _poly);
                _c ^= _d;
                
                _GET_X_PAK(i)[ii] = _POLY_GET_LOWER_FEC_INT(_c);
                _GET_X_PAK(i+1)[ii] = _POLY_GET_UPPER_FEC_INT(_c);
            }
            
            if (i == num_y_missing - 1) {
                _poly_t _c = _POLY_INT_TO_POLY(tmp_recovered_ints[i]);
                _poly_t _d = _c >> 16;
                _d = _POLY_SINGLE_MUL(_d, _poly);
                _c ^= _d;
                _d >>= 16;
                _d = _POLY_SINGLE_MUL(_d, _poly);
                _c ^= _d;

                _GET_X_PAK(i)[ii] = _POLY_TO_FEC_INT(_c);
            }
        }
#elif defined(FEC_HAS_128_INT_VEC) || defined(FEC_HAS_64_INT_VEC)
        {
            for (i = 0; i < num_y_missing; i++) {
                tmp_recovered_ints[i] = (u16x16){ones_pak_ii, 0};
            }

            for (j = 0; j < num_y_present + num_x_present; j++) {
                fec_int_t xy_j = present_y[j];
                fec_int_t xy_pak_arr_j_val = _GET_XY_PAK(j)[ii];
#if defined(FEC_HAS_128_INT_VEC) || (defined(FEC_HAS_64_INT_VEC) && !(defined(__x86_64__) || defined(__i386__)))
                __fec_rx_col_op_5(tmp_recovered_ints, missing_y, num_y_missing, inv_arr, xy_pak_arr_j_val, xy_j);
#elif (defined(__x86_64__) || defined(__i386__)) && defined(__SSE__)
                fec_rx_col_op_sse((__m128*)tmp_recovered_ints, missing_y, num_y_missing, inv_arr, xy_pak_arr_j_val, xy_j);
#elif (defined(__x86_64__) || defined(__i386__)) && defined(__MMX__)
                fec_rx_col_op_mmx((__m64*)tmp_recovered_ints, missing_y, num_y_missing, inv_arr, xy_pak_arr_j_val, xy_j);
#endif          
            }

#if defined(FEC_HAS_128_INT_VEC)
            u32x8 shifts1 = {0, 2, 4, 6, 8, 10, 12, 14};
            u32x8 shifts2 = {1, 3, 5, 7, 9, 11, 13, 15};

            for (i = 0; i < num_y_missing; i++) {
                u16x16 res = tmp_recovered_ints[i];
                u32x8 res1 = (((u32x8)res) & 0xFFFF) << shifts1;
                u32x8 res2 = (((u32x8)res) >> 16) << shifts2;

                u32x8 res3 = res1 ^ res2;
                u32x4 res4 = (u32x4)(my_mm256_extracti128_si256((__m256i)res3, 0) ^ my_mm256_extracti128_si256((__m256i)res3, 1));
                u32x4 res5 = res4 ^ my_mm_shuffle_epi32(res4, 1, 2, 3, 0);
                u32x4 res6 = res5 ^ my_mm_shuffle_epi32(res5, 2, 3, 0, 1); // res6 holds in each cell the int32 res

                u32x4 res7 = (res6 >> 16) << ((u32x4){0, 1, 3, 5}); // mul upper by POLY_G
                
                u32x4 res8 = res7 ^ my_mm_shuffle_epi32(res7, 1, 2, 3, 0);
                u32x4 res9 = res8 ^ my_mm_shuffle_epi32(res8, 2, 3, 0, 1);

                u32x4 res10 = (res9 >> 16) << ((u32x4){0, 1, 3, 5}); // mul upper by POLY_G
                
                u32x4 res11 = res10 ^ my_mm_shuffle_epi32(res10, 1, 2, 3, 0);
                u32x4 res12 = res11 ^ my_mm_shuffle_epi32(res11, 2, 3, 0, 1);

                _GET_X_PAK(i)[ii] = ((u16x8)(res6 ^ res9 ^ res12))[0];

            }
#else

            for (i = 0; i < num_y_missing; i++) {

                uint32_t res = 0;
                int jj;
                for (jj = 0; jj < 16; jj++) {
                    res ^= ((uint32_t)(((uint16_t*)tmp_recovered_ints)[i*16 + jj])) << jj;
                }

                uint16_t res4 = res;
                uint32_t carry = res >> 16;

                carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

                res4 ^= carry;

                carry = carry >> 16;
                
                res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

                _GET_X_PAK(i)[ii] = res4;
            }

#if defined(__x86_64__) || defined(__i386__)
            _mm_empty();
#endif
#endif
        }
#elif defined(FEC_HAS_64BIT)
    {
        for (i = 0; i < num_y_missing; i++) {
            ((uint64_t*)&tmp_recovered_ints[i])[0] = 0;
            ((uint64_t*)&tmp_recovered_ints[i])[1] = 0;
            ((uint64_t*)&tmp_recovered_ints[i])[2] = 0;
            ((uint64_t*)&tmp_recovered_ints[i])[3] = (((uint64_t)ones_pak_ii) << (64-16));
        }
        
        for (j = 0; j < num_y_present + num_x_present; j++) {
            fec_int_t xy_j = present_y[j];
            fec_int_t xy_pak_arr_j_val = _GET_XY_PAK(j)[ii];
            fec_rx_col_op_reg64((uint64_t*)tmp_recovered_ints, missing_y, num_y_missing, inv_arr, xy_pak_arr_j_val, xy_j);
        }

        for (i = 0; i < num_y_missing; i++) {
            uint32_t res = 0;
            int jj;
            for (jj = 0; jj < 16; jj++) {
                res ^= ((uint32_t)(((uint16_t*)tmp_recovered_ints)[i*16 + jj])) << (16 - 1 - jj);
            }

            uint16_t res4 = res;
            uint32_t carry = res >> 16;

            carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

            res4 ^= carry;

            carry = carry >> 16;
            
            res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

            _GET_X_PAK(i)[ii] = res4;
        }
    }

#elif defined(FEC_HAS_32BIT)
    {
        for (i = 0; i < num_y_missing; i++) {
            unsigned int jj;
            ((uint32_t*)&tmp_recovered_ints[i])[0] = ones_pak_ii;
            for (jj = 1; jj < 8; jj++) {
                ((uint32_t*)&tmp_recovered_ints[i])[jj] = 0;
            }
        }
        
        for (j = 0; j < num_y_present + num_x_present; j++) {
            fec_int_t xy_j = present_y[j];
            fec_int_t xy_pak_arr_j_val = _GET_XY_PAK(j)[ii];
            fec_rx_col_op_reg32((uint32_t*)tmp_recovered_ints, missing_y, num_y_missing, inv_arr, xy_pak_arr_j_val, xy_j);
        }

        for (i = 0; i < num_y_missing; i++) {
            uint32_t res = 0;
            int jj;
            for (jj = 0; jj < 16; jj++) {
                res ^= ((uint32_t)(((uint16_t*)tmp_recovered_ints)[i*16 + jj])) << jj;
            }

            uint16_t res4 = res;
            uint32_t carry = res >> 16;

            carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

            res4 ^= carry;

            carry = carry >> 16;
            
            res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

            _GET_X_PAK(i)[ii] = res4;
        }
    }

#else
        for (i = 0; i < num_y_missing; i++) {
            tmp_recovered_ints[i] = ones_pak_ii;
        }

        for (j = 0; j < num_y_present + num_x_present; j++) {
            fec_int_t xy_j = present_y[j];
            fec_int_t xy_pak_arr_j_ii = _GET_XY_PAK(j)[ii];
            for (i = 0; i < num_y_missing; i++) {
                fec_int_t missing_y_i = missing_y[i];
                tmp_recovered_ints[i] ^= poly_mul(xy_pak_arr_j_ii, _fec_inv(inv_cache, poly_add(xy_j, missing_y_i)));
            }
        }
        for (i = 0; i < num_y_missing; i++) {
            _GET_X_PAK(i)[ii] = ((fec_int_t*)tmp_recovered_ints)[i];
        }
#endif
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.4--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t pi_yx_div_yy_i = 1;
        y_i = missing_y[i];
        for (j = 0; j < num_x_present; j++) {
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, poly_add(y_i, present_x[j]));
        }
        for (j = 0; j < num_y_missing; j++) {
            if(j == i) {
                continue;
            }
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, _fec_inv(inv_cache, poly_add(y_i, missing_y[j])));
        }

        unaligend_fec_int_t *pak = _GET_X_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_yx_div_yy_i);
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.5--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    for (i = 0; i < num_y_present; i++) {
        y_i = present_y[i];

        fec_int_t inv_pi_ycomp_y_div_ycomp_x_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i, _fec_inv(inv_cache, poly_add(y_i, missing_y[j])));
        }
        for (j = 0; j < num_x_present; j++) {
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i,  poly_add(y_i, present_x[j]));
        }

        unaligend_fec_int_t* pak = _GET_Y_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], inv_pi_ycomp_y_div_ycomp_x_i);
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.6--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    // {
    //     uint16_t aaa = rand();
    //     uint16_t bbb = rand();
    //     for(i = 0; i < n*pak_len*num_y_missing; i++) {
    //         aaa = poly_mul(aaa, bbb);
    //         bbb = aaa ^ bbb;
    //     }
    //     printf("---%d---\n", aaa);
    // }

    for (i = 0; i < num_y_missing; i++) {
        present_x[i] = missing_y[i];
    }

    //blabla(num_y_missing, n, pak_len, has_one_row);

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.7--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

reorder_packets:
#ifdef PERF_DEBUG
    start_time = get_timestamp();
#endif
    for (i = 0; i < n; i++) {
        while (rx_state->pak_xy_arr[i] != i) {
            fec_int_t idx = rx_state->pak_xy_arr[i];

            {
                fec_int_t tmp;
                tmp = rx_state->pak_xy_arr[idx];
                rx_state->pak_xy_arr[idx] = idx;
                rx_state->pak_xy_arr[i] = tmp;
            }

#ifndef FEC_USER_GIVEN_BUFFER
            {
                unaligend_fec_int_t *tmp;
                tmp = rx_state->pak_arr[idx];
                rx_state->pak_arr[idx] = rx_state->pak_arr[i];
                rx_state->pak_arr[i] = tmp;
            }
#else
            memswap(&rx_state->pak_buffer[i*pak_len], &rx_state->pak_buffer[idx*pak_len], pak_len);
#endif
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.8--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    return FEC_STATUS_SUCCESS;
}

#endif

#endif


fec_status_t fec_rx_get_needed_inv_cache_size(fec_rx_state_t *rx_state, fec_idx_t* n_k_1) {
    if (!_fec_can_recover(rx_state)) {
        return FEC_STATUS_MORE_PACKETS_NEEDED;
    }

    if (rx_state->num_info != rx_state->n && rx_state->max_x != 0) {
        // we return the minimal "size" needed to pass to fec_inv_cache_init_raw
        // we do this so the user can compare it to what he initialized with
        *n_k_1 = (1 << (log2_ceil(((fec_idx_t)rx_state->max_x) + 1) - 1)) + 1;
    } else {
        *n_k_1 = 0;
    }

    return FEC_STATUS_SUCCESS;
}


#ifndef FEC_USER_GIVEN_BUFFER
void** fec_rx_get_info_paks(const fec_rx_state_t *rx_state) {
#ifndef FEC_LARGE_K
    return (void**)rx_state->info_paks;
#else
    return (void**)rx_state->pak_arr;
#endif
}
#endif

