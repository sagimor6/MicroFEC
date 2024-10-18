
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#ifdef PERF_DEBUG
#include <stdio.h>
#ifdef _WIN32
#include <realtimeapiset.h>
#include <processthreadsapi.h>
#include <profileapi.h>
#else
#include <time.h>
#endif

#if defined(__x86_64__) && !defined(__SSE2__)
#define PRINT_TS_DIFF(str) do { \
    end_time = get_timestamp(); \
    uint64_t __diff = end_time - start_time; \
    printf(str " %lu.%09lu\n", __diff / 1000000000, __diff % 1000000000); \
    start_time = get_timestamp(); \
} while(0)
#else
#define PRINT_TS_DIFF(str) do { \
    end_time = get_timestamp(); \
    uint64_t __diff = end_time - start_time; \
    printf(str " %f\n", __diff / ((double)1000000000)); \
    start_time = get_timestamp(); \
} while(0)
#endif
#endif

#include "micro_fec.h"
#include "fec_common.h"

#ifdef FEC_DO_ENDIAN_SWAP
#include <byteswap.h>
#endif

#ifdef PERF_DEBUG
static uint64_t get_timestamp() {
#ifdef _WIN32
    // uint64_t t;
    // QueryThreadCycleTime(GetCurrentThread(), &t);
    // t /= 4;
    // return t;

    uint64_t t, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&t);
    return (uint64_t)(((double)t * 1000000000) / freq);
#else
    struct timespec tp = {0};
    //CHECK(clock_gettime(CLOCK_MONOTONIC, &tp) == 0);
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
    return (tp.tv_sec*1000000000ULL) + tp.tv_nsec;
#endif
}
#endif


#ifdef FEC_DO_ENDIAN_SWAP
#define fec_bswap(x) (bswap_16(x))
#else
#define fec_bswap(x) (x)
#endif

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

#if defined(FEC_HAS_CLMUL32)
#include "opt_templates/clmul.h"
#endif

static fec_int_t PERF_DEBUG_ATTRS poly_mul(fec_int_t a, fec_int_t b) {
#if defined(FEC_HAS_CLMUL32)

    _poly_t _a = _POLY_1VAL(a);
    _poly_t _b = _POLY_1VAL(b);
    _poly_t _poly = _POLY_1VAL(POLY_G);

    _poly_t _c;
    _c = _POLY_CLMUL(_a, _b);
    _poly_t _d;
    _d = _c >> 16;
    _d = _POLY_CLMUL(_d, _poly);
    _c ^= _d;
    _d >>= 16;
    _d = _POLY_CLMUL(_d, _poly);
    _c ^= _d;
    
    return _POLY_EXTRACT(_c, uint16_t, 0);
#elif defined(FEC_HAS_32BIT)
    size_t i;
    uint32_t res = 0;

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        res ^= (((uint32_t)a) << i) & ((uint32_t)(int32_t)(((int16_t)(b << (15 - i))) >> 15));
    }

    uint16_t res2 = res;
    uint32_t carry = res >> 16;

    carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

    res2 ^= carry;

    carry = carry >> 16;
    
    res2 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

    return res2;
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

#ifdef FEC_DO_ENDIAN_SWAP
static void PERF_DEBUG_ATTRS mem_bswap(unaligned_fec_int_t* arr, size_t size) {
    // TODO: i am hoping compiler will optimize this

    size_t i;
    for (i = 0; i < size; i++) {
        fec_int_t val = arr[i];
        arr[i] = bswap_16(val);
    }
}
#endif

// TODO: decided not to copy in reverse endian at user given buffer at fec_rx_add_pak

// #ifdef FEC_USER_GIVEN_BUFFER
// static void PERF_DEBUG_ATTRS fec_memcpy_bswap(unaligned_fec_int_t* dst, const unaligned_fec_int_t* src, size_t size) {
// #ifdef FEC_DO_ENDIAN_SWAP
//     size_t i;
//     for (i = 0; i < size; i++) {
//         dst[i] = bswap_16(src[i]);
//     }
// #else
//     memcpy(dst, src, size * sizeof(fec_int_t));
// #endif
// }
// #endif

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
    size_t tmp_pak_elems = pak_len;
#ifdef PERF_TX_BLOCK_SIZE
    tmp_pak_elems = MIN(tmp_pak_elems, PERF_TX_BLOCK_SIZE);
#endif
    tx_state->tmp_pak = malloc(tmp_pak_elems * sizeof(tx_state->tmp_pak[0]) + CALC_EXTRA_FOR_ALIGNED_MALLOC(tx_state->tmp_pak[0]));
    if (tx_state->tmp_pak == NULL) {
        fec_tx_destroy(tx_state);
        return FEC_STATUS_OUT_OF_MEMORY;
    }
#endif

    return FEC_STATUS_SUCCESS;
}

void fec_tx_reset(fec_tx_state_t *tx_state) {
    memset(tx_state->paks, 0, tx_state->n * sizeof(tx_state->paks[0]));
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

#ifndef FEC_USER_GIVEN_BUFFER
    MALLOC_ATTR(pak_arr, n);
#else
    rx_state->pak_buffer = (unaligned_fec_int_t*)dest_buf;
#endif

    CALLOC_ATTR(received_paks_bitmap, (n + k + (8-1))/8);
    MALLOC_ATTR(pak_xy_arr, n);

    MALLOC_ATTR(missing_y, MIN(k, n));

#ifdef FEC_MIN_MEM
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT) && defined(PERF_RX_BLOCK_SIZE)
    MALLOC_ATTR(tmp_recovered_block, MIN(MIN(k, n), PERF_RX_BLOCK_SIZE));
#endif
    MALLOC_ATTR(tmp_recovered_ints, MIN(k, n));
#else
    MALLOC_ATTR(pak_multiplier, n);
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT)
    size_t tmp_pak_elems = pak_len;
#ifdef PERF_RX_BLOCK_SIZE
    tmp_pak_elems = MIN(tmp_pak_elems, PERF_RX_BLOCK_SIZE);
#endif
    MALLOC_ATTR(tmp_pak, tmp_pak_elems);
#endif
#endif

#undef MALLOC_ATTR
#undef CALLOC_ATTR

    rx_state->num_info = 0;
    rx_state->num_redundant = 0;
    rx_state->max_x = n - 1;
#ifndef FEC_USER_GIVEN_BUFFER
    rx_state->ones_pak = NULL;
#else
    rx_state->has_one_pak = false;
#endif

    return FEC_STATUS_SUCCESS;
}

void fec_rx_reset(fec_rx_state_t *rx_state) {
    if (rx_state->n == 0) {
        return;
    }
    memset(rx_state->received_paks_bitmap, 0, ((rx_state->n + rx_state->k + (8-1))/8) * sizeof(rx_state->received_paks_bitmap[0]));
    rx_state->num_info = 0;
    rx_state->num_redundant = 0;
    rx_state->max_x = rx_state->n - 1;
#ifndef FEC_USER_GIVEN_BUFFER
    rx_state->ones_pak = NULL;
#else
    rx_state->has_one_pak = false;
#endif
}

void fec_rx_destroy(fec_rx_state_t *rx_state) {
#define FREE_ATTR(name) \
    if (rx_state->name != NULL) { \
        free(rx_state->name); \
    }
    
#ifndef FEC_USER_GIVEN_BUFFER
    FREE_ATTR(pak_arr);
#endif
    FREE_ATTR(received_paks_bitmap);
    FREE_ATTR(pak_xy_arr);

    FREE_ATTR(missing_y);

#ifdef FEC_MIN_MEM
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT) && defined(PERF_RX_BLOCK_SIZE)
    FREE_ATTR(tmp_recovered_block);
#endif
    FREE_ATTR(tmp_recovered_ints);
#else
    FREE_ATTR(pak_multiplier);
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT)
    FREE_ATTR(tmp_pak);
#endif
#endif
#undef FREE_ATTR
}

fec_status_t fec_tx_add_info_pak(fec_tx_state_t *tx_state, const void* pak, fec_idx_t idx) {
    if (idx >= tx_state->n) {
        return FEC_STATUS_INVALID_PARAMS;
    }
    tx_state->paks[idx] = (const unaligned_fec_int_t*)pak;
    return FEC_STATUS_SUCCESS;
}

static bool _fec_can_recover(const fec_rx_state_t *rx_state) {
#if !defined(FEC_USER_GIVEN_BUFFER)
    return (rx_state->num_info + rx_state->num_redundant + (rx_state->ones_pak != NULL) >= rx_state->n);
#else
    return (rx_state->num_info + rx_state->num_redundant >= rx_state->n); // we count the ones pak in the redundant
#endif
}


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
        rx_state->pak_arr[rx_state->num_info] = (unaligned_fec_int_t*)pak;
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
        rx_state->pak_arr[n - 1 - rx_state->num_redundant] = (unaligned_fec_int_t*)pak;
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

#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_TX_OPT)

static void fec_tx_init_perf_arr(fec_perf_int_t* restrict out_pak, size_t pak_len) {
    memset(out_pak, 0, pak_len*sizeof(out_pak[0]));
}

#define LEN_PARAM_TYPE size_t
#define INPUT_ARGS const unaligned_fec_int_t* restrict pak
#if defined(FEC_DO_ENDIAN_SWAP) && defined(FEC_HAS_CLMUL32)
#define READ_INPUT(j) bswap_16(pak[j])
#else
#define READ_INPUT(j) pak[j]
#endif
#define OUTPUT_ARGS unaligned_fec_int_t* restrict out_pak
#define WRITE_OUTPUT(j, val) out_pak[j] = fec_bswap(val)
#define INIT_FUNC_NAME -1
#define FMA_FUNC_NAME fec_tx_col_op
#define NORM_FUNC_NAME fec_tx_col_perf_to_norm
#include "opt_templates/perf.c"
#undef FMA_FUNC_NAME
#undef NORM_FUNC_NAME
#undef INIT_FUNC_NAME
#undef LEN_PARAM_TYPE
#undef INPUT_ARGS
#undef READ_INPUT
#undef OUTPUT_ARGS
#undef WRITE_OUTPUT
#endif

// TODO: idx can be fec_int_t
fec_status_t fec_tx_get_redundancy_pak(const fec_tx_state_t *tx_state, const fec_inv_cache_t *inv_cache, fec_idx_t idx, void *pak) {
    fec_idx_t n = tx_state->n;
    size_t pak_len = tx_state->pak_len;
    fec_idx_t i;
    size_t j;
    unaligned_fec_int_t* out_pak = (unaligned_fec_int_t*)pak;
    const unaligned_fec_int_t* restrict * restrict paks = tx_state->paks;
    
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_TX_OPT)
    fec_perf_int_t* restrict tmp_pak = tx_state->tmp_pak;
    tmp_pak = (fec_perf_int_t*)(((uintptr_t)tmp_pak + __alignof__(fec_perf_int_t) - 1) & (-__alignof__(fec_perf_int_t)));
    const fec_int_t* restrict inv_arr = inv_cache->inv_arr - 1;
#endif

    if (n == 0 || pak_len == 0) {
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

    if (idx == 0) {
        memset(out_pak, 0, pak_len*sizeof(fec_int_t));

        for (i = 0; i < n; i++) {
            const unaligned_fec_int_t* pak = paks[i];
            for (j = 0; j < pak_len; j++) {
                out_pak[j] = poly_add(out_pak[j], pak[j]);
                // we don't need endian swap here because bswap(bswap(a)^bswap(b)) = a^b
            }
        }
        return FEC_STATUS_SUCCESS;
    }

#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_TX_OPT)
#ifdef PERF_TX_BLOCK_SIZE
    for (j = 0; j < pak_len; j += PERF_TX_BLOCK_SIZE) {
        size_t cur_pak_len = MIN(PERF_TX_BLOCK_SIZE, pak_len - j);
#else
    {
        j = 0;
        size_t cur_pak_len = pak_len;
#endif
        fec_tx_init_perf_arr(tmp_pak, cur_pak_len);
        for (i = 0; i < n; i++) {
            fec_int_t a_i = inv_arr[poly_add(n + idx - 1, i)];
            const unaligned_fec_int_t* pak = paks[i];
            fec_tx_col_op(tmp_pak, cur_pak_len, a_i, &pak[j]);
        }
#if defined(FEC_DO_ENDIAN_SWAP) && !defined(FEC_HAS_CLMUL32)
        // before normalizing, bswap(a)*b = bswap(a*b), bswap(a^b) = bswap(a)^bswap(b)
        mem_bswap((fec_int_t*)tmp_pak, (sizeof(tmp_pak[0])/(sizeof(fec_int_t))) * cur_pak_len);
#endif
        fec_tx_col_perf_to_norm(tmp_pak, cur_pak_len, &out_pak[j]);
    }

#if !defined(FEC_HAS_128_INT_VEC) && defined(FEC_HAS_64_INT_VEC) && (defined(__x86_64__) || defined(__i386__))
    _mm_empty();
#endif

#else
    memset(out_pak, 0, pak_len*sizeof(fec_int_t));

    for (i = 0; i < n; i++) {
        fec_int_t a_i = _fec_inv(inv_cache, poly_add(n + idx - 1, i));
        const unaligned_fec_int_t* pak = paks[i];
        for (j = 0; j < pak_len; j++) {
            // if (idx == 0) {
            //     out_pak[j] = poly_add(out_pak[j], pak[j]);
            // } else {
                out_pak[j] = poly_add(out_pak[j], poly_mul(fec_bswap(pak[j]), a_i));
            // }
        }
    }
#if defined(FEC_DO_ENDIAN_SWAP)
    mem_bswap(out_pak, pak_len);
#endif
#endif

    return FEC_STATUS_SUCCESS;
}

#ifndef FEC_USER_GIVEN_BUFFER
#define _GET_X_PAK(idx) (x_pak_arr[(idx)])
#define _GET_Y_PAK(idx) (y_pak_arr[(idx)])
#define _GET_XY_PAK(idx) _GET_Y_PAK(idx)
#define _NORM_FMA_PARAMS_DEC unaligned_fec_int_t* restrict* restrict x_pak_arr
#define _NORM_FMA_PARAMS x_pak_arr
#else
#define _GET_X_PAK(idx) (&x_paks_buf[(idx)*pak_len])
#define _GET_Y_PAK(idx) (&y_paks_buf[(idx)*pak_len])
#define _GET_XY_PAK(idx) _GET_Y_PAK(idx)
#define _NORM_FMA_PARAMS_DEC unaligned_fec_int_t* restrict x_paks_buf, size_t pak_len
#define _NORM_FMA_PARAMS x_paks_buf, pak_len
#endif

#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT)
#ifdef FEC_MIN_MEM
#define LEN_PARAM_TYPE fec_idx_t
#define INIT_FUNC_ARGS fec_int_t val
#define INIT_FUNC_READ_INPUT(j) val
#define INPUT_ARGS const fec_int_t* restrict missing_y, const fec_int_t* restrict inv_arr, fec_int_t pak_xy
#define READ_INPUT(j) inv_arr[poly_add(pak_xy, missing_y[j])]
#ifdef PERF_RX_BLOCK_SIZE
#define OUTPUT_ARGS fec_int_t* restrict tmp_recovered_ints
#define WRITE_OUTPUT(j, val) tmp_recovered_ints[j] = val
#else
#define OUTPUT_ARGS size_t ii, _NORM_FMA_PARAMS_DEC
#define WRITE_OUTPUT(j, val) _GET_X_PAK(j)[ii] = val
#endif
#define INIT_FUNC_NAME fec_rx_col_init
#define FMA_FUNC_NAME fec_rx_col_op
#define NORM_FUNC_NAME fec_rx_col_perf_to_norm
#include "opt_templates/perf.c"
#undef FMA_FUNC_NAME
#undef NORM_FUNC_NAME
#undef INIT_FUNC_NAME
#undef LEN_PARAM_TYPE
#undef INPUT_ARGS
#undef READ_INPUT
#undef OUTPUT_ARGS
#undef WRITE_OUTPUT
#undef INIT_FUNC_ARGS
#undef INIT_FUNC_READ_INPUT
#else
#define LEN_PARAM_TYPE size_t
#define INIT_FUNC_ARGS const fec_int_t* restrict ones_pak
#define INIT_FUNC_READ_INPUT(j) ones_pak[j]
#define INPUT_ARGS const fec_int_t* restrict pak
#define READ_INPUT(j) pak[j]
#define OUTPUT_ARGS unaligned_fec_int_t* restrict out_pak
#define WRITE_OUTPUT(j, val) out_pak[j] = val
#define INIT_FUNC_NAME fec_rx_col_init2
#define FMA_FUNC_NAME fec_rx_col_op2
#define NORM_FUNC_NAME fec_rx_col_perf_to_norm2
#include "opt_templates/perf.c"
#undef FMA_FUNC_NAME
#undef NORM_FUNC_NAME
#undef INIT_FUNC_NAME
#undef LEN_PARAM_TYPE
#undef INPUT_ARGS
#undef READ_INPUT
#undef OUTPUT_ARGS
#undef WRITE_OUTPUT
#undef INIT_FUNC_ARGS
#undef INIT_FUNC_READ_INPUT
#endif
#endif

#ifdef FEC_USER_GIVEN_BUFFER
static void PERF_DEBUG_ATTRS memswap(unaligned_fec_int_t* a, unaligned_fec_int_t* b, size_t size) {
    // TODO: i am hoping compiler will optimize this

    size_t i;
    for (i = 0; i < size; i++) {
        fec_int_t tmp = a[i];
        a[i] = b[i];
        b[i] = tmp;
    }
}
#endif

fec_status_t fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state, const fec_inv_cache_t *inv_cache) {
    fec_idx_t n = rx_state->n;
    size_t pak_len = rx_state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present;
    fec_int_t *present_x;
    fec_int_t *present_y;
    fec_int_t *missing_y;
    unaligned_fec_int_t *ones_pak = NULL;
    fec_idx_t num_y_present = rx_state->num_info;
#ifndef FEC_MIN_MEM
    fec_int_t *pak_multiplier = rx_state->pak_multiplier;
#endif

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
#ifdef FEC_DO_ENDIAN_SWAP
    for (i = 0; i < n; i++) {
        mem_bswap(rx_state->pak_arr[i], pak_len);
    }
#endif
#else
    has_one_row = rx_state->has_one_pak;
    num_x_present = num_y_missing - has_one_row;
    if (has_one_row && rx_state->ones_pak_idx != 0) {
        // we need to move the ones pak
        memswap(&rx_state->pak_buffer[(n - 1 - rx_state->ones_pak_idx)*pak_len], &rx_state->pak_buffer[(n - 1)*pak_len], pak_len);
        rx_state->pak_xy_arr[n - 1 - rx_state->ones_pak_idx] = rx_state->pak_xy_arr[n - 1];
        //rx_state->ones_pak_idx = 0;
    }
#ifdef FEC_DO_ENDIAN_SWAP
    mem_bswap(rx_state->pak_buffer, n * pak_len);
#endif
#endif

    present_x = &rx_state->pak_xy_arr[num_y_present];
    present_y = rx_state->pak_xy_arr;
    missing_y = rx_state->missing_y;

#ifndef FEC_USER_GIVEN_BUFFER
    unaligned_fec_int_t** x_pak_arr = &rx_state->pak_arr[num_y_present];
    unaligned_fec_int_t** y_pak_arr = rx_state->pak_arr;
    ones_pak = rx_state->ones_pak;
#else
    unaligned_fec_int_t* x_paks_buf = &rx_state->pak_buffer[num_y_present * pak_len];
    unaligned_fec_int_t* y_paks_buf = rx_state->pak_buffer;
    ones_pak = &rx_state->pak_buffer[(n - 1)*pak_len];
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
    PRINT_TS_DIFF("---1.1---");
#endif

#ifdef FEC_MIN_MEM
#define _Y_X_MUL_START_IDX 0
#else
#define _Y_X_MUL_START_IDX 1
#endif
    
    for (i = 0; i < num_y_present; i++) {
        y_i = present_y[i];

        fec_int_t pi_ycomp_y_div_ycomp_x_i = 1;
        for (j = _Y_X_MUL_START_IDX; j < num_y_missing; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, poly_add(y_i, missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, _fec_inv(inv_cache, poly_add(y_i, present_x[j])));
        }

#ifndef FEC_MIN_MEM
        pak_multiplier[i] = pi_ycomp_y_div_ycomp_x_i;
#else
        unaligned_fec_int_t* pak = _GET_Y_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_ycomp_y_div_ycomp_x_i);
        }
#endif
    }

#ifdef PERF_DEBUG
    PRINT_TS_DIFF("---1.2---");
#endif

    for (i = 0; i < num_x_present; i++) {
        fec_int_t pi_xy_div_xx_i = 1;
        fec_int_t x_i = present_x[i];
        for (j = _Y_X_MUL_START_IDX; j < num_y_missing; j++) {
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, poly_add(x_i, missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, _fec_inv(inv_cache, poly_add(x_i, present_x[j])));
        }

#ifndef FEC_MIN_MEM
        pak_multiplier[num_y_present + i] = pi_xy_div_xx_i;
#else
        unaligned_fec_int_t *pak = _GET_X_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_xy_div_xx_i);
        }
#endif
    }

#ifdef PERF_DEBUG
    PRINT_TS_DIFF("---1.3---");
#endif

#if !defined(FEC_MIN_MEM) || (!defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT))
    const fec_int_t* inv_arr = inv_cache->inv_arr - 1;
#endif

#ifndef FEC_MIN_MEM

#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT)
    fec_perf_int_t *tmp_pak = (fec_perf_int_t*)(((uintptr_t)rx_state->tmp_pak + __alignof__(fec_perf_int_t) - 1) & (-__alignof__(fec_perf_int_t)));
#endif

    fec_int_t y_j;
    for (j = 0; j < num_y_missing; j++) {
        y_j = missing_y[j];

#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT)
        unaligned_fec_int_t* rec_pak = _GET_X_PAK(j);
#ifdef PERF_RX_BLOCK_SIZE
        for (ii = 0; ii < pak_len; ii += PERF_RX_BLOCK_SIZE) {
            size_t cur_block_size = MIN(PERF_RX_BLOCK_SIZE, pak_len - ii);
#else
        {
            ii = 0;
            size_t cur_block_size = pak_len;
#endif
            
            if (has_one_row) {
                fec_rx_col_init2(tmp_pak, cur_block_size, ones_pak);
            } else {
                memset(tmp_pak, 0, sizeof(tmp_pak[0])*cur_block_size);
            }

            for (i = 0; i < n - has_one_row; i++) {
                fec_rx_col_op2(tmp_pak, cur_block_size, pak_multiplier[i], _GET_XY_PAK(i));
            }
            
            fec_rx_col_perf_to_norm2(tmp_pak, cur_block_size, &rec_pak[ii]);
        }
#else
        unaligned_fec_int_t* rec_pak = _GET_X_PAK(j);
        if (j != num_y_missing - 1 || !has_one_row) {
            fec_int_t mult = pak_multiplier[num_y_present];
            for (ii = 0; ii < pak_len; ii++) {
                rec_pak[ii] = poly_mul(rec_pak[ii], mult);
            }

            for(i = 0; i < n - has_one_row; i++) {
                if (i == num_y_present) {
                    continue;
                }
                unaligned_fec_int_t* pak = _GET_XY_PAK(i);
                fec_int_t mult = pak_multiplier[i];
                for(ii = 0; ii < pak_len; ii++) {
                    rec_pak[ii] ^= poly_mul(pak[ii], mult);
                }
            }
            
            if (has_one_row) {
                for(ii = 0; ii < pak_len; ii++) {
                    rec_pak[ii] ^= ones_pak[ii];
                }
            }
        } else {
            for (i = 0; i < n - 1; i++) {
                unaligned_fec_int_t* pak = _GET_XY_PAK(i);
                fec_int_t mult = pak_multiplier[i];
                for(ii = 0; ii < pak_len; ii++) {
                    rec_pak[ii] ^= poly_mul(pak[ii], mult);
                }
            }
        }
#endif

        // fec_int_t pi_yx_div_yy = 1;
        // for (i = j; i < num_x_present; i++) {
        //     pi_yx_div_yy = poly_mul(pi_yx_div_yy, poly_add(y_j, present_x[i]));
        // }
        // for (i = j + 1; i < num_y_missing; i++) {
        //     pi_yx_div_yy = poly_mul(pi_yx_div_yy, _fec_inv(inv_cache, poly_add(y_j, missing_y[i])));
        // }

        // fec_int_t pi_ycomp_y_div_ycomp_x = 1;
        // for (i = j + 1; i < num_y_missing; i++) {
        //     pi_ycomp_y_div_ycomp_x = poly_mul(pi_ycomp_y_div_ycomp_x, poly_add(y_j, missing_y[i]));
        // }
        // for (i = j + 1; i < num_x_present; i++) {
        //     pi_ycomp_y_div_ycomp_x = poly_mul(pi_ycomp_y_div_ycomp_x, _fec_inv(inv_cache, poly_add(y_j, present_x[i])));
        // }

        if (j != num_y_missing - 1) {
            fec_int_t x_j = present_x[0];
            fec_int_t y_j_1 = missing_y[j + 1];

            for (i = 0; i < num_y_present; i++) {
                fec_int_t y_i = present_y[i];
                fec_int_t mul_correct = inv_arr[poly_add(y_i, y_j_1)];
                mul_correct = poly_mul(mul_correct, poly_add(y_i, x_j));
                pak_multiplier[i] = poly_mul(pak_multiplier[i], mul_correct);
            }

            pak_multiplier[num_y_present] = poly_mul(poly_add(y_j, x_j), inv_arr[poly_add(y_j, y_j_1)]);

            for (i = 1; i < num_x_present; i++) {
                fec_int_t x_i = present_x[i];
                fec_int_t mul_correct = inv_arr[poly_add(x_i, y_j_1)];
                mul_correct = poly_mul(mul_correct, poly_add(x_i, x_j));
                pak_multiplier[num_y_present + i] = poly_mul(pak_multiplier[num_y_present + i], mul_correct);
            }

            present_x[0] = y_j;
            present_x++;
            num_x_present--;
            num_y_present++;
        } else {
            if (!has_one_row) {
                fec_int_t x_j = present_x[0];
                for (i = 0; i < num_y_present; i++) {
                    fec_int_t y_i = present_y[i];
                    pak_multiplier[i] = poly_mul(pak_multiplier[i], poly_add(y_i, x_j));
                }

                pak_multiplier[num_y_present] = poly_add(y_j, x_j);

                present_x[0] = y_j;
                present_x++;
                num_y_present++;
                num_x_present--;
            } else {
                present_x[0] = y_j;
                num_y_present++;
            }
        }
        
        // for (ii = 0; ii < pak_len; ii++) {
        //     rec_pak[ii] = poly_mul(rec_pak[ii], pi_yx_div_yy);
        // }
    }
#else

#if defined(_FEC_NO_OPT) || defined(_FEC_NO_RX_OPT)
    fec_int_t* tmp_recovered_ints = rx_state->tmp_recovered_ints;
#elif defined(PERF_RX_BLOCK_SIZE)
    fec_perf_int_t* tmp_recovered_block = (fec_perf_int_t*)(((uintptr_t)rx_state->tmp_recovered_block + __alignof__(fec_perf_int_t) - 1) & (-__alignof__(fec_perf_int_t)));
    fec_int_t* tmp_recovered_ints = rx_state->tmp_recovered_ints;
#else
    fec_perf_int_t* tmp_recovered_ints = (fec_perf_int_t*)(((uintptr_t)rx_state->tmp_recovered_ints + __alignof__(fec_perf_int_t) - 1) & (-__alignof__(fec_perf_int_t)));
#endif

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

        //     res = __fec_rx_row_op((const unaligned_fec_int_t* const*)y_pak_arr, num_y_present + num_x_present, ii, present_y, missing_y_i, inv_cache->inv_arr, ones_pak_ii);

        //     tmp_recovered_ints[i] = res;
        // }

        // for (i = 0; i < num_y_missing; i++) {
        //     _GET_X_PAK(i)[ii] = ((fec_int_t*)tmp_recovered_ints)[i];
        // }


#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_RX_OPT)
#ifndef PERF_RX_BLOCK_SIZE
        fec_rx_col_init(tmp_recovered_ints, num_y_missing, ones_pak_ii);

        for (j = 0; j < num_y_present + num_x_present; j++) {
            fec_int_t xy_j = present_y[j];
            fec_int_t xy_pak_arr_j_val = _GET_XY_PAK(j)[ii];
            fec_rx_col_op(tmp_recovered_ints, num_y_missing, xy_pak_arr_j_val, missing_y, inv_arr, xy_j);
        }

        fec_rx_col_perf_to_norm(tmp_recovered_ints, num_y_missing, ii, _NORM_FMA_PARAMS);
#else
        for (i = 0; i < num_y_missing; i += PERF_RX_BLOCK_SIZE) {
            fec_idx_t cur_block_size = MIN(PERF_RX_BLOCK_SIZE, num_y_missing - i);
            fec_rx_col_init(tmp_recovered_block, cur_block_size, ones_pak_ii);

            for (j = 0; j < num_y_present + num_x_present; j++) {
                fec_int_t xy_j = present_y[j];
                fec_int_t xy_pak_arr_j_val = _GET_XY_PAK(j)[ii];
                fec_rx_col_op(tmp_recovered_block, cur_block_size, xy_pak_arr_j_val, &missing_y[i], inv_arr, xy_j);
            }

            fec_rx_col_perf_to_norm(tmp_recovered_block, cur_block_size, &tmp_recovered_ints[i]);
        }
        for (i = 0; i < num_y_missing; i++) {
            _GET_X_PAK(i)[ii] = tmp_recovered_ints[i];
        }
#endif
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
#endif

#if !defined(FEC_HAS_128_INT_VEC) && defined(FEC_HAS_64_INT_VEC) && (defined(__x86_64__) || defined(__i386__))
    _mm_empty();
#endif

#ifdef PERF_DEBUG
    PRINT_TS_DIFF("---1.4---");
#endif

#ifndef FEC_MIN_MEM
    for (i = n - num_y_missing; i < n - has_one_row; i++) {
        unaligned_fec_int_t *pak = _GET_XY_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pak_multiplier[i]);
        }
    }
#else
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

        unaligned_fec_int_t *pak = _GET_X_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_yx_div_yy_i);
        }
    }

#ifdef PERF_DEBUG
    PRINT_TS_DIFF("---1.5---");
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

        unaligned_fec_int_t* pak = _GET_Y_PAK(i);
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], inv_pi_ycomp_y_div_ycomp_x_i);
        }
    }

#ifdef PERF_DEBUG
    PRINT_TS_DIFF("---1.6---");
#endif

    for (i = 0; i < num_y_missing; i++) {
        present_x[i] = missing_y[i];
    }
#endif

#ifdef PERF_DEBUG
    PRINT_TS_DIFF("---1.7---");
#endif

#ifdef FEC_DO_ENDIAN_SWAP
#ifndef FEC_USER_GIVEN_BUFFER
    for (i = 0; i < n; i++) {
        mem_bswap(y_pak_arr[i], pak_len);
    }
#else
    mem_bswap(y_paks_buf, n * pak_len);
#endif
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
                unaligned_fec_int_t *tmp;
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
    PRINT_TS_DIFF("---1.8---");
#endif

    return FEC_STATUS_SUCCESS;
}


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
    return (void**)rx_state->pak_arr;
}
#endif

