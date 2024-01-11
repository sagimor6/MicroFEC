
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

//#define PERF_DEBUG 
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

typedef uint8_t u8x8 __attribute__ ((vector_size (8)));
typedef uint16_t u16x4 __attribute__ ((vector_size (8)));
typedef uint32_t  u32x2 __attribute__ ((vector_size (8)));
typedef uint64_t  u64x1 __attribute__ ((vector_size (8)));

#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#endif

#if (defined(__PCLMUL__) && defined(__SSE2__)) || defined(_FEC_NO_OPT) || !(defined(__x86_64__) || defined(__i386__))
#define _FEC_USE_POLY_MUL
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

static fec_int_t PERF_DEBUG_ATTRS poly_mul(fec_int_t a, fec_int_t b) {
#if defined(__PCLMUL__) && defined(__SSE2__)
    // uint32_t res;

    // asm(
    //     "movd %k[_a], %%xmm0\n"
    //     "movd %k[_b], %%xmm1\n"
    //     "movd %k[_poly], %%xmm2\n"
    //     "PCLMULLQLQDQ %%xmm1, %%xmm0\n"
    //     "movq %%xmm0, %%xmm1\n"
    //     "PSRLD $%c[_shift], %%xmm1\n"
    //     "PCLMULLQLQDQ %%xmm2, %%xmm1\n"
    //     "XORPS %%xmm1, %%xmm0\n"
    //     "PSRLD $%c[_shift], %%xmm1\n"
    //     "PCLMULLQLQDQ %%xmm2, %%xmm1\n"
    //     "XORPS %%xmm1, %%xmm0\n"
    //     "movd %%xmm0, %k[_out]\n"
    // : [_out] "=r" (res)
    // : [_a] "r" ((uint32_t)a), [_b] "r" ((uint32_t)b), [_poly] "r" ((uint32_t)POLY_G), [_shift] "i" (sizeof(fec_int_t)*8)
    // : "%xmm0", "%xmm1", "%xmm2"
    // );
    // return (fec_int_t)res;


    // v128 _a;
    // _a.u64[0] = a;
    // v128 _b;
    // _b.u64[0] = b;
    // v128 _poly;
    // _poly.u64[0] = POLY_G;

    // v128 _c;
    // _c.mm = _mm_clmulepi64_si128(_a.mm, _b.mm, 0);
    // v128 _d;
    // _d.u32 = (_c.u32 >> 16);
    // _d.mm = _mm_clmulepi64_si128(_d.mm, _poly.mm, 0);
    // _c.u16 ^= _d.u16;
    // _d.u32 = _d.u32 >> 16;
    // _d.mm = _mm_clmulepi64_si128(_d.mm, _poly.mm, 0);
    // _c.u16 ^= _d.u16;
    // return _c.u16[0];

    u64x2 _a;
    _a[0] = a;
    u64x2 _b;
    _b[0] = b;
    u64x2 _poly;
    _poly[0] = POLY_G;

    u16x8 _c;
    _c = (u16x8)_mm_clmulepi64_si128((__m128i)_a, (__m128i)_b, 0);
    u32x4 _d;
    _d = (((u32x4)_c) >> 16);
    _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    _c ^= (u16x8)_d;
    _d >>= 16;
    _d = (u32x4)_mm_clmulepi64_si128((__m128i)_d, (__m128i)_poly, 0);
    _c ^= (u16x8)_d;
    return _c[0];

#elif defined(__ARM_FEATURE_AES)
    // TODO: FEAT_PMULL is needed in processor
#ifdef __aarch64__
#define _poly_t poly64_t
#else
#define _poly_t uint32_t
#endif
    _poly_t _c;
    _c = (_poly_t)vmull_p64(a, b);
    _poly_t _d;
    _d = (_c >> 16);
    _d = (_poly_t)vmull_p64(_d, POLY_G);
    _c ^= _d;
    _d >>= 16;
    _d = (_poly_t)vmull_p64(_d, POLY_G);
    _c ^= _d;
    return (fec_int_t)_c;
#undef _poly_t
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

static inline fec_int_t _fec_inv(const fec_state_t *state, fec_int_t a) {
    return state->inv_arr[a - 1];
}

bool fec_init(fec_state_t *state, fec_idx_t n, fec_idx_t k, size_t pak_len) {
    state->n = n;
    state->k = k;
    state->pak_len = pak_len;

    // TODO: i can deal with this case
    if (n == 0 || k == 0) {
        return false;
    }

    fec_idx_t n_k_1 = n + k - 1;
    if (n_k_1 < n || n_k_1 > (((fec_idx_t)1)<<(sizeof(fec_int_t)*8))) {
        return false;
    }

    if (n_k_1 != 0) {
        size_t arr_size = (1 << log2_ceil(n_k_1)) - 1;
        state->inv_arr = malloc(sizeof(state->inv_arr[0]) * arr_size);
        if (state->inv_arr == NULL) {
            return false;
        }

        fec_int_t i;
        for (i = 0; i < arr_size; i++) {
            state->inv_arr[i] = poly_inv(i + 1);
        }

    } else {
        state->inv_arr = NULL;
    }

    return true;
}

void fec_destroy(fec_state_t *state) {
    if (state->inv_arr != NULL) {
        free(state->inv_arr);
    }
}

bool fec_tx_init(fec_tx_state_t *tx_state, fec_state_t *state) {
    tx_state->state = state;
    tx_state->paks = calloc(state->n, sizeof(tx_state->paks[0]));
    if (tx_state->paks == NULL) {
        return false;
    }

    return true;
}

void fec_tx_destroy(fec_tx_state_t *tx_state) {
    if (tx_state->paks != NULL) {
        free(tx_state->paks);
    }
}

bool fec_rx_init(fec_rx_state_t *rx_state, fec_state_t *state) {
    memset(rx_state, 0, sizeof(*rx_state));

#define MALLOC_ATTR(name, size) \
    rx_state->name = malloc(sizeof(rx_state->name[0]) * (size)); \
    if (rx_state->name == NULL) { \
        fec_rx_destroy(rx_state); \
        return false; \
    }
#define CALLOC_ATTR(name, size) \
    rx_state->name = calloc((size), sizeof(rx_state->name[0])); \
    if (rx_state->name == NULL) { \
        fec_rx_destroy(rx_state); \
        return false; \
    }

    fec_idx_t n = state->n; 
    fec_idx_t k = state->k;

    rx_state->state = state;
#ifndef FEC_LARGE_K
    CALLOC_ATTR(info_paks, n);
    CALLOC_ATTR(redundancy_paks, k);
    MALLOC_ATTR(present_x, MIN(k - 1, n));
#else
    CALLOC_ATTR(received_paks_bitmap, (n + k + (8-1))/8);
    MALLOC_ATTR(pak_arr, n);
    MALLOC_ATTR(pak_xy_arr, n);
#endif

    MALLOC_ATTR(missing_y, MIN(k, n));

#ifdef FEC_MIN_MEM
    MALLOC_ATTR(tmp_recovered_ints, ALIGN_UP(MIN(k, n), _FEC_ALIGN_SIZE_VAL/sizeof(fec_int_t)));
#else
    MALLOC_ATTR(pi_xy_div_xx, MIN(k - 1, n));
    MALLOC_ATTR(pi_yx_div_yy, MIN(k, n));
    MALLOC_ATTR(tmp_vec_redundancy, MIN(k - 1, n));
    MALLOC_ATTR(tmp_vec_info, state->n - 1);
    MALLOC_ATTR(present_y, state->n - 1);
    MALLOC_ATTR(pi_ycomp_y_div_ycomp_x, state->n - 1);
#endif

#undef MALLOC_ATTR
#undef CALLOC_ATTR

    rx_state->num_info = 0;
    rx_state->num_redundant = 0;
#ifdef FEC_LARGE_K
    rx_state->ones_pak = NULL;
#endif
    return true;
}

void fec_rx_reset(fec_rx_state_t *rx_state) {
#ifndef FEC_LARGE_K
    memset(rx_state->info_paks, 0, rx_state->state->n * sizeof(rx_state->info_paks[0]));
    memset(rx_state->redundancy_paks, 0, rx_state->state->k * sizeof(rx_state->redundancy_paks[0]));
#else
    memset(rx_state->received_paks_bitmap, 0, ((rx_state->state->n + rx_state->state->k + (8-1))/8) * sizeof(rx_state->received_paks_bitmap[0]));
#endif
    rx_state->num_info = 0;
    rx_state->num_redundant = 0;
#ifdef FEC_LARGE_K
    rx_state->ones_pak = NULL;
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
    FREE_ATTR(received_paks_bitmap);
    FREE_ATTR(pak_arr);
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

bool fec_tx_add_info_pak(fec_tx_state_t *tx_state, const void* pak, fec_idx_t idx) {
    if (idx >= tx_state->state->n) {
        return false;
    }
    tx_state->paks[idx] = (const unaligend_fec_int_t*)pak;
    return true;
}

#ifndef FEC_LARGE_K

bool fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_idx_t idx, bool *can_recover, bool *discard_pak) {
    fec_int_t n = rx_state->state->n;
    bool _can_recover;
    bool _discard_pak = true;

    if (idx >= n + rx_state->state->k) {
        return false;
    }

    _can_recover = (rx_state->num_info + rx_state->num_redundant >= n);

    if (_can_recover) {
        goto valid_pak_end;
    }

    if (idx < n) {
        if (rx_state->info_paks[idx] != NULL) {
            goto valid_pak_end;
        }
        rx_state->info_paks[idx] = (unaligend_fec_int_t*)pak;
        rx_state->num_info++;
    } else {
        if (rx_state->redundancy_paks[idx - n] != NULL) {
            goto valid_pak_end;
        }
        rx_state->redundancy_paks[idx - n] = (unaligend_fec_int_t*)pak;
        if (idx != n) {
            bool has_one_row = (rx_state->redundancy_paks[0] != NULL);
            rx_state->present_x[rx_state->num_redundant - has_one_row] = idx - 1;
        }
        rx_state->num_redundant++;
    }

    _discard_pak = false;
    _can_recover = (rx_state->num_info + rx_state->num_redundant >= n);
    
valid_pak_end:
    if (can_recover) {
        *can_recover = _can_recover;
    }
    if (discard_pak) {
        *discard_pak = _discard_pak;
    }
    return true;
}

#else

bool fec_rx_is_pak_needed(fec_rx_state_t *rx_state, fec_idx_t idx, bool *can_recover, bool *discard_pak) {
    fec_int_t n = rx_state->state->n;
    bool _can_recover;
    bool _discard_pak;

    if (idx >= n + rx_state->state->k) {
        return false;
    }

    _can_recover = (rx_state->num_info + rx_state->num_redundant + (rx_state->ones_pak != NULL) >= n);

    // if can recover or duplicate
    _discard_pak = _can_recover || (rx_state->received_paks_bitmap[idx / 8] & (1<<(idx & (8-1)))) != 0;

    *can_recover = _can_recover;
    *discard_pak = _discard_pak;
    return true;
}

bool fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_idx_t idx, bool *can_recover, bool *discard_pak) {
    fec_int_t n = rx_state->state->n;
    bool _can_recover;
    bool _discard_pak = true;

    if (idx >= n + rx_state->state->k) {
        return false;
    }

    _can_recover = (rx_state->num_info + rx_state->num_redundant + (rx_state->ones_pak != NULL) >= n);

    if (_can_recover) {
        goto valid_pak_end;
    }

    if ((rx_state->received_paks_bitmap[idx / 8] & (1<<(idx & (8-1)))) != 0) {
        goto valid_pak_end; // already received
    }

    _discard_pak = false;

    rx_state->received_paks_bitmap[idx / 8] |= (1<<(idx & (8-1)));

    if (idx < n) {
        rx_state->pak_arr[rx_state->num_info] = (unaligend_fec_int_t*)pak;
        rx_state->pak_xy_arr[rx_state->num_info] = idx;
        rx_state->num_info++;
    } else if (idx == n) {
        rx_state->ones_pak = pak;
    } else {
        rx_state->pak_arr[n - 1 - rx_state->num_redundant] = (unaligend_fec_int_t*)pak;
        rx_state->pak_xy_arr[n - 1 - rx_state->num_redundant] = idx - 1;
        rx_state->num_redundant++;
    }

    _can_recover = (rx_state->num_info + rx_state->num_redundant + (rx_state->ones_pak != NULL) >= n);

valid_pak_end:
    if (can_recover) {
        *can_recover = _can_recover;
    }
    if (discard_pak) {
        *discard_pak = _discard_pak;
    }
    return true;
}

#endif

bool fec_tx_get_redundancy_pak(const fec_tx_state_t *tx_state, fec_idx_t idx, void *pak) {
    const fec_state_t *state = tx_state->state;
    fec_idx_t n = state->n;
    size_t pak_len = state->pak_len;
    fec_idx_t i;
    size_t j;
    unaligend_fec_int_t* out_pak = (unaligend_fec_int_t*)pak;

    if (idx >= tx_state->state->k) {
        return false;
    }

    // TODO: this does not affect the time complexity of the function, but maybe remember number in struct
    for (i = 0; i < n; i++) {
        if (tx_state->paks[i] == NULL) {
            return false;
        }
    }

    // for (j = 0; j < pak_len; j++) {
    //     fec_int_t res = 0;
    //     for (i = 0; i < n; i++) {
    //         if(idx == 0) {
    //             res = poly_add(res, tx_state->paks[i][j]);
    //         } else {
    //             fec_int_t a_i = poly_add(n + idx - 1, i);
    //             res = poly_add(res, poly_mul(tx_state->paks[i][j], _fec_inv(state, a_i)));
    //         }
    //     }
    //     out_pak[j] = res;
    // }

    memset(out_pak, 0, pak_len*sizeof(fec_int_t));

    if (idx == 0) {
        for (i = 0; i < n; i++) {
            const unaligend_fec_int_t* pak = tx_state->paks[i];
            for (j = 0; j < pak_len; j++) {
                out_pak[j] = poly_add(out_pak[j], pak[j]);
            }
        }
        return true;
    }

    for (i = 0; i < n; i++) {
        fec_int_t a_i = _fec_inv(state, poly_add(n + idx - 1, i));
        const unaligend_fec_int_t* pak = tx_state->paks[i];
        for (j = 0; j < pak_len; j++) {
            // if (idx == 0) {
            //     out_pak[j] = poly_add(out_pak[j], pak[j]);
            // } else {
                out_pak[j] = poly_add(out_pak[j], poly_mul(pak[j], a_i));
            // }
        }
    }

    return true;
}

#ifndef FEC_MIN_MEM
#ifndef FEC_LARGE_K
bool fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state) {
    const fec_state_t *state = rx_state->state;
    fec_idx_t n = state->n;
    size_t pak_len = state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present = 0;

    fec_idx_t i, j;
    fec_idx_t y_i, y_j;
    size_t ii;

    if (rx_state->num_redundant < num_y_missing) {
        return false;
    }

    if (num_y_missing == 0) {
        return true;
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
            res = poly_mul(res, _fec_inv(state, poly_add(rx_state->present_x[i], rx_state->present_x[j])));
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
            res = poly_mul(res, _fec_inv(state, poly_add(rx_state->missing_y[i], rx_state->missing_y[j])));
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
            res = poly_mul(res, _fec_inv(state, poly_add(y_i, rx_state->present_x[j])));
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
                res = poly_add(res, poly_mul(rx_state->tmp_vec_info[j], _fec_inv(state, poly_add(y_j, rx_state->missing_y[i]))));
                j++;
            }
            if (has_one_row) {
                res = poly_add(res, tmp_vec_1s);
            }
            //for (x_j = 0, j = 0; x_j < k - 1; x_j++) {
            for(j = 0; j < num_x_present; j++) {
                res = poly_add(res, poly_mul(rx_state->tmp_vec_redundancy[j], _fec_inv(state, poly_add(rx_state->present_x[j], rx_state->missing_y[i]))));
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

    return true;
}

#else

bool fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state) {
    const fec_state_t *state = rx_state->state;
    fec_idx_t n = state->n;
    size_t pak_len = state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present = 0;
    fec_int_t *present_x;

    fec_idx_t i, j;
    fec_idx_t y_i, y_j;
    size_t ii;

    if (rx_state->num_info + rx_state->num_redundant + (rx_state->ones_pak != NULL) < n) {
        return false;
    }

    if (num_y_missing == 0) {
        return true;
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
            res = poly_mul(res, _fec_inv(state, poly_add(present_x[i], present_x[j])));
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
            res = poly_mul(res, _fec_inv(state, poly_add(rx_state->missing_y[i], rx_state->missing_y[j])));
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
            res = poly_mul(res, _fec_inv(state, poly_add(y_i, present_x[j])));
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
                res = poly_add(res, poly_mul(rx_state->tmp_vec_info[j], _fec_inv(state, poly_add(y_j, rx_state->missing_y[i]))));
            }
            if (has_one_row) {
                res = poly_add(res, tmp_vec_1s);
            }
            for(j = 0; j < num_x_present; j++) {
                res = poly_add(res, poly_mul(rx_state->tmp_vec_redundancy[j], _fec_inv(state, poly_add(present_x[j], rx_state->missing_y[i]))));
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

    return true;
}

#endif

#else

#ifndef FEC_LARGE_K

bool fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state) {
    const fec_state_t *state = rx_state->state;
    fec_idx_t n = state->n;
    size_t pak_len = state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present = 0;

    fec_idx_t i, j;
    fec_idx_t y_i, y_j;
    size_t ii;

    if (rx_state->num_redundant < num_y_missing) {
        return false;
    }

    if (num_y_missing == 0) {
        return true;
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
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, _fec_inv(state, poly_add(y_i, rx_state->present_x[j])));
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
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, _fec_inv(state, poly_add(rx_state->present_x[i], rx_state->present_x[j])));
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
                res = poly_add(res, poly_mul(rx_state->info_paks[y_j][ii], _fec_inv(state, poly_add(y_j, rx_state->missing_y[i]))));
                j++;
            }
            if (has_one_row) {
                res = poly_add(res, rx_state->redundancy_paks[0][ii]);
            }
            //for (x_j = 0, j = 0; x_j < k - 1; x_j++) {
            for(j = 0; j < num_x_present; j++) {
                res = poly_add(res, poly_mul(rx_state->redundancy_paks[rx_state->present_x[j] - n + 1][ii], _fec_inv(state, poly_add(rx_state->present_x[j], rx_state->missing_y[i]))));
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
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, _fec_inv(state, poly_add(rx_state->missing_y[i], rx_state->missing_y[j])));
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
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i, _fec_inv(state, poly_add(y_i, rx_state->missing_y[j])));
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

    return true;
}

#else


static void PERF_DEBUG_ATTRS __fec_rx_col_op(fec_int_t* recovered, const fec_int_t *missing_y, fec_idx_t num_y_missing, const fec_int_t* inv_arr, fec_int_t pak_val, fec_int_t pak_xy) {

    inv_arr = inv_arr - 1;

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
        u16x16 cc = poly_mul16(aa, bb);

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
    
}


#ifdef PERF_DEBUG
static uint64_t get_timestamp() {
    struct timespec tp = {0};
    //CHECK(clock_gettime(CLOCK_MONOTONIC, &tp) == 0);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp);
    return (tp.tv_sec*1000000000ULL) + tp.tv_nsec;
}
#endif

bool fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state) {
    const fec_state_t *state = rx_state->state;
    fec_idx_t n = state->n;
    size_t pak_len = state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row = (rx_state->ones_pak != NULL);
    fec_idx_t num_x_present;
    fec_int_t *present_x;
    fec_int_t *present_y;
    fec_int_t *missing_y;
    unaligend_fec_int_t** x_pak_arr;
    unaligend_fec_int_t** y_pak_arr;
    fec_idx_t num_y_present = rx_state->num_info;

    fec_idx_t i, j;
    fec_idx_t y_i, y_j;
    size_t ii;

    if (rx_state->num_info + rx_state->num_redundant + has_one_row < n) {
        return false;
    }

    if (num_y_missing == 0) {
        goto reorder_packets;
    }

    num_x_present = num_y_missing - has_one_row;

    present_x = &rx_state->pak_xy_arr[n - num_x_present];
    present_y = rx_state->pak_xy_arr;
    missing_y = rx_state->missing_y;
    x_pak_arr = &rx_state->pak_arr[n - num_x_present];
    y_pak_arr = rx_state->pak_arr;

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
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, _fec_inv(state, poly_add(y_i, present_x[j])));
        }
        for (ii = 0; ii < pak_len; ii++) {
            y_pak_arr[i][ii] = poly_mul(y_pak_arr[i][ii], pi_ycomp_y_div_ycomp_x_i);
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.2--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    for (i = 0; i < num_x_present; i++) {
        fec_int_t pi_xy_div_xx_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, poly_add(present_x[i], missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, _fec_inv(state, poly_add(present_x[i], present_x[j])));
        }

        unaligend_fec_int_t *pak = x_pak_arr[i];
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_xy_div_xx_i);
        }
    }

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.3--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    fec_int_t* tmp_recovered_ints = rx_state->tmp_recovered_ints;

    for (ii = 0; ii < pak_len; ii++) {
        

        fec_int_t ones_pak_ii = 0;
        if (has_one_row) {
            ones_pak_ii = rx_state->ones_pak[ii];
        }

        // for (i = 0; i < num_y_missing; i++) {

        //     fec_int_t res;

        //     res = ones_pak_ii;

        //     fec_int_t missing_y_i = missing_y[i];

        //     for (j = 0; j < num_y_present; j++) {
        //         y_j = present_y[j];
        //         res = poly_add(res, poly_mul(y_pak_arr[j][ii], _fec_inv(state, poly_add(y_j, missing_y_i))));
        //     }
        //     for(j = 0; j < num_x_present; j++) {
        //         res = poly_add(res, poly_mul(x_pak_arr[j][ii], _fec_inv(state, poly_add(present_x[j], missing_y_i))));
        //     }

        //     tmp_recovered_ints[i] = res;
        // }

        for (i = 0; i < num_y_missing; i++) {
            tmp_recovered_ints[i] = ones_pak_ii;
        }            

        for (j = 0; j < num_y_present; j++) {
            y_j = present_y[j];
            fec_int_t y_pak_arr_j_ii = y_pak_arr[j][ii];
            // for (i = 0; i < num_y_missing; i++) {
            //     fec_int_t missing_y_i = missing_y[i];
            //     tmp_recovered_ints[i] ^= poly_mul(y_pak_arr_j_ii, _fec_inv(state, poly_add(y_j, missing_y_i)));
            // }
            __fec_rx_col_op(tmp_recovered_ints, missing_y, num_y_missing, state->inv_arr, y_pak_arr_j_ii, y_j);
        }
        for(j = 0; j < num_x_present; j++) {
            fec_int_t x_j = present_x[j];
            fec_int_t x_pak_arr_j_ii = x_pak_arr[j][ii];
            // for (i = 0; i < num_y_missing; i++) {
            //     fec_int_t missing_y_i = missing_y[i];
            //     tmp_recovered_ints[i] ^= poly_mul(x_pak_arr_j_ii, _fec_inv(state, poly_add(x_j, missing_y_i)));
            // }
            __fec_rx_col_op(tmp_recovered_ints, missing_y, num_y_missing, state->inv_arr, x_pak_arr_j_ii, x_j);
        }

        for (i = 0; i < num_y_missing; i++) {
            if (i == 0 && has_one_row) {
                rx_state->ones_pak[ii] = tmp_recovered_ints[i];
            } else {
                x_pak_arr[i - has_one_row][ii] = tmp_recovered_ints[i];
            }
        }
    }

#ifdef _FEC_USE_POLY_MUL4_MMX
    _mm_empty();
#endif

#ifdef PERF_DEBUG
    end_time = get_timestamp();
    printf("---1.4--- %f\n", (end_time - start_time)/((double)1000000000));
    start_time = get_timestamp();
#endif

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t pi_yx_div_yy_i = 1;
        for (j = 0; j < num_x_present; j++) {
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, poly_add(missing_y[i], present_x[j]));
        }
        for (j = 0; j < num_y_missing; j++) {
            if(j == i) {
                continue;
            }
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, _fec_inv(state, poly_add(missing_y[i], missing_y[j])));
        }

        for (ii = 0; ii < pak_len; ii++) {
            if (i == 0 && has_one_row) {
                rx_state->ones_pak[ii] = poly_mul(rx_state->ones_pak[ii], pi_yx_div_yy_i);
            } else {
                x_pak_arr[i - has_one_row][ii] = poly_mul(x_pak_arr[i - has_one_row][ii], pi_yx_div_yy_i);
            }
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
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i, _fec_inv(state, poly_add(y_i, missing_y[j])));
        }
        for (j = 0; j < num_x_present; j++) {
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i,  poly_add(y_i, present_x[j]));
        }
        for (ii = 0; ii < pak_len; ii++) {
            y_pak_arr[i][ii] = poly_mul(y_pak_arr[i][ii], inv_pi_ycomp_y_div_ycomp_x_i);
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

    if (has_one_row) {
        rx_state->pak_xy_arr[n - num_x_present - 1] = rx_state->missing_y[0];
        rx_state->pak_arr[n - num_x_present - 1] = rx_state->ones_pak;
    }
    for (i = 0; i < num_x_present; i++) {
        rx_state->pak_xy_arr[n - num_x_present + i] = rx_state->missing_y[i + has_one_row];
    }

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

    return true;
}

#endif

#endif


void** fec_rx_get_info_paks(const fec_rx_state_t *rx_state) {
#ifndef FEC_LARGE_K
    return (void**)rx_state->info_paks;
#else
    return (void**)rx_state->pak_arr;
#endif
}

