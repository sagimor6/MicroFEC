#ifndef __MICRO_FEC_H__
#define __MICRO_FEC_H__

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

typedef uint16_t fec_int_t;
typedef uint32_t fec_idx_t;

#ifdef FEC_MIN_MEM
#undef FEC_MIN_MEM
#endif
#ifdef FEC_LARGE_K
#undef FEC_LARGE_K
#endif

#define FEC_LARGE_K
#define FEC_MIN_MEM

typedef fec_int_t __attribute__((aligned(1))) unaligend_fec_int_t;

typedef enum {
    FEC_STATUS_SUCCESS,
    FEC_STATUS_OUT_OF_MEMORY,
    FEC_STATUS_INVALID_PARAMS,
    FEC_STATUS_INV_CACHE_TOO_SMALL,
    FEC_STATUS_MORE_PACKETS_NEEDED,
    FEC_STATUS_CAN_DROP_DUP_PAK,
    FEC_STATUS_CAN_DROP_ALREADY_RECOVERABLE,
} fec_status_t;

typedef struct {
    // n+k-1 <= (1<<(sizeof(fec_int_t)*8))

    // 0,...,n+k-2 are the numbers
    // 1,...,1<<ceil(log2(n+k-1))

    fec_idx_t inv_arr_sz; // ceil(log2(n+k-1))
    fec_int_t* inv_arr; // upper power of two closest to n+k-2

} fec_inv_cache_t;


#if !defined(_FEC_NO_OPT)

#if ((defined(__x86_64__) || defined(__i386__)) && defined(__PCLMUL__)) || \
    ((defined(__aarch64__) || defined(__arm__)) && defined(__ARM_FEATURE_AES)) || /* TODO: FEAT_PMULL is needed in processor */ \
    (defined(__sparc__) && defined(__VIS) && __VIS >= 0x300) || \
    (defined(__riscv__) && (defined(__riscv_zbc) || defined(__riscv_zbkc)) && __riscv_xlen == 64) || \
    0
#define FEC_HAS_CLMUL64
#endif

#if defined(FEC_HAS_CLMUL64) || \
    (defined(__riscv__) && (defined(__riscv_zbc) || defined(__riscv_zbkc)) && __riscv_xlen == 32) || \
    0
#define FEC_HAS_CLMUL32
#endif

#if ((defined(__x86_64__) || defined(__i386__)) && defined(__SSE2__)) || \
    ((defined(__aarch64__) || defined(__arm__)) && __ARM_NEON) || \
    (defined(__mips__) && __mips_msa) || \
    (defined(__powerpc__) && __ALTIVEC__) || \
    (defined(__riscv__) && __riscv_vector && __riscv_v_min_vlen >= 128/8 && __riscv_v_elen >= 32) /* TODO: check this */ || \
    0
#define FEC_HAS_128_INT_VEC
#endif

#if defined(FEC_HAS_128_INT_VEC) || \
    ((defined(__x86_64__) || defined(__i386__)) && defined(__MMX__)) || \
    (defined(__sparc__) && defined(__VIS) && __VIS >= 0x100) || \
    (defined(__riscv__) && __riscv_vector && __riscv_v_min_vlen >= 64/8 && __riscv_v_elen >= 32) /* TODO: check this */ || \
    0
#define FEC_HAS_64_INT_VEC
#endif

#if defined(__LP64__) || \
    defined(__x86_64__) || \
    defined(__aarch64__) || \
    defined(__mips64) || \
    defined(__powerpc64__) || \
    (defined(__riscv__) && __riscv_xlen == 64) || \
    0
#define FEC_HAS_64BIT
#endif

#if __INT_MAX__ > __SHRT_MAX__
#define FEC_HAS_32BIT
#endif

// very weird case because x86_64 imples SSE2
// but if x86_64 with only mmx, then we prefer regular 64bit impl
#if defined(__x86_64__) && !defined(FEC_HAS_128_INT_VEC) && defined(FEC_HAS_64_INT_VEC) && !defined(__SSE__)
#undef FEC_HAS_64_INT_VEC
#endif

// #if defined(FEC_HAS_CLMUL64)
// #define FEC_USE_CLMUL64
// #elif 

// #endif

#endif

// #undef FEC_HAS_CLMUL64
// #undef FEC_HAS_CLMUL32
// #undef FEC_HAS_128_INT_VEC
// #undef FEC_HAS_64_INT_VEC
// #undef FEC_HAS_64BIT
// #undef FEC_HAS_32BIT

#if defined(FEC_HAS_CLMUL32)
typedef uint32_t fec_perf_int_t;
#elif defined(FEC_HAS_64_INT_VEC)
#if defined(FEC_HAS_128_INT_VEC) || !(defined(__x86_64__) || defined(__i386__))
typedef uint16_t __attribute__ ((vector_size (32))) fec_perf_int_t;
#elif (defined(__x86_64__) || defined(__i386__)) && defined(__SSE__)
typedef uint16_t __attribute__ ((vector_size (32))) __attribute__((aligned(16))) fec_perf_int_t;
#elif (defined(__x86_64__) || defined(__i386__)) && defined(__MMX__)
typedef uint16_t __attribute__ ((vector_size (32))) __attribute__((aligned(8))) fec_perf_int_t;
#endif
#elif defined(FEC_HAS_64BIT)
typedef uint64_t fec_perf_int_t[4];
#elif defined(FEC_HAS_32BIT)
typedef uint32_t fec_perf_int_t[8];
#else
typedef uint16_t fec_perf_int_t;
#endif

typedef struct {
    fec_idx_t n;
    size_t pak_len;

    const unaligend_fec_int_t** paks; // size = n
#if !defined(_FEC_NO_OPT) && !defined(_FEC_NO_TX_OPT)
    fec_perf_int_t* tmp_pak;
#endif
} fec_tx_state_t;

//#define FEC_USER_GIVEN_BUFFER

// typedef uint16_t __attribute__ ((vector_size (32))) fec_2int_t;
// typedef uint32_t fec_2int_t;
// typedef uint16_t fec_2int_t;

typedef struct {
    fec_idx_t n;
    fec_idx_t k;
    size_t pak_len;

    fec_int_t max_x; // maximum k - 1

#ifndef FEC_LARGE_K
    unaligend_fec_int_t** info_paks; // size = n
    unaligend_fec_int_t** redundancy_paks; // size = real k
    fec_int_t *present_x; // size = k - 1
#else
    uint8_t* received_paks_bitmap; // size = (n + real k)/8
#ifndef FEC_USER_GIVEN_BUFFER
    unaligend_fec_int_t** pak_arr; // size = n
    unaligend_fec_int_t* ones_pak;
#else
    unaligend_fec_int_t *pak_buffer; // given by user (size = n*pak_len)
    fec_int_t ones_pak_idx;
    bool has_one_pak; // TODO: actually x_i cannot be 0, cause n > 0, so ones_pak_idx can point to elem in pak_xy_arr with 0.
    // TODO: we also have this info from the bitfield
#endif
    fec_int_t *pak_xy_arr; // size = n
#endif
    fec_idx_t num_info;
    fec_idx_t num_redundant;

    fec_int_t *missing_y; // size = k

#ifdef FEC_MIN_MEM
    // min mem uses this:
    fec_perf_int_t *tmp_recovered_ints; // size = k
#else
    // regular uses this:
    fec_int_t *pi_xy_div_xx; // size = k - 1
    fec_int_t *pi_yx_div_yy; // size = k

    fec_int_t *present_y; // size = n - 1

    fec_int_t *pi_ycomp_y_div_ycomp_x; // size = n - 1

    fec_int_t *tmp_vec_info; // size = n - 1
    fec_int_t *tmp_vec_redundancy; // size = k - 1
#endif
} fec_rx_state_t;

#ifdef _FEC_DO_EXPORTS
#define EXPORT __attribute__((visibility("default")))
#else
#define EXPORT 
#endif

EXPORT fec_status_t fec_inv_cache_init(fec_inv_cache_t *inv_cache, fec_idx_t n, fec_idx_t k);
EXPORT fec_status_t fec_inv_cache_init_raw(fec_inv_cache_t *inv_cache, fec_idx_t n_k_1);
EXPORT void fec_inv_cache_destroy(fec_inv_cache_t *inv_cache);

EXPORT fec_status_t fec_tx_init(fec_tx_state_t *tx_state, fec_idx_t n, size_t pak_len);
EXPORT fec_status_t fec_tx_add_info_pak(fec_tx_state_t *tx_state, const void* pak, fec_idx_t idx);
EXPORT fec_status_t fec_tx_get_redundancy_pak(const fec_tx_state_t *tx_state, const fec_inv_cache_t *inv_cache, fec_idx_t idx, void *pak);
EXPORT void fec_tx_destroy(fec_tx_state_t *tx_state);

#ifndef FEC_USER_GIVEN_BUFFER
EXPORT fec_status_t fec_rx_init(fec_rx_state_t *rx_state, fec_idx_t n, fec_idx_t k, size_t pak_len);
#else
EXPORT fec_status_t fec_rx_init(fec_rx_state_t *rx_state, fec_idx_t n, fec_idx_t k, size_t pak_len, void* dest_buf);
#endif
EXPORT fec_status_t fec_rx_is_pak_needed(fec_rx_state_t *rx_state, fec_idx_t idx);
EXPORT fec_status_t fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_idx_t idx);
EXPORT fec_status_t fec_rx_get_needed_inv_cache_size(fec_rx_state_t *rx_state, fec_idx_t* n_k_1);
EXPORT fec_status_t fec_rx_fill_missing_paks(const fec_rx_state_t *rx_state, const fec_inv_cache_t *inv_cache);
#ifndef FEC_USER_GIVEN_BUFFER
EXPORT void** fec_rx_get_info_paks(const fec_rx_state_t *rx_state);
#endif
EXPORT void fec_rx_reset(fec_rx_state_t *rx_state);
EXPORT void fec_rx_destroy(fec_rx_state_t *rx_state);

#undef EXPORT

#endif // __MICRO_FEC_H__
