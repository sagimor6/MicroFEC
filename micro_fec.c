
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "micro_fec.h"

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

static fec_int_t poly_mul(fec_int_t a, fec_int_t b) {
#if defined(__PCLMUL__) && defined(__SSE2__)
    fec_int_t res;

    asm(
        ".intel_syntax noprefix\n"
        "movd xmm0, %k[_a]\n"
        "movd xmm1, %k[_b]\n"
        "movd xmm2, %k[_poly]\n"
        "PCLMULLQLQDQ xmm0, xmm1\n"
        "movq xmm1, xmm0\n"
        "PSRLD xmm1, %c[_shift]\n"
        "PCLMULLQLQDQ xmm1, xmm2\n"
        "XORPS xmm0, xmm1\n"
        "PSRLD xmm1, %c[_shift]\n"
        "PCLMULLQLQDQ xmm1, xmm2\n"
        "XORPS xmm0, xmm1\n"
        "movd %k[_out], xmm0\n"
        ".att_syntax prefix\n"
    : [_out] "=r" (res)
    : [_a] "r" (a), [_b] "r" (b), [_poly] "r" (POLY_G), [_shift] "i" (sizeof(fec_int_t)*8)
    : "xmm0", "xmm1", "xmm2"
    );
    return res;
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
    MALLOC_ATTR(tmp_recovered_ints, MIN(k, n));
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

    for (j = 0; j < pak_len; j++) {
        fec_int_t res = 0;
        for (i = 0; i < n; i++) {
            if(idx == 0) {
                res = poly_add(res, tx_state->paks[i][j]);
            } else {
                fec_int_t a_i = poly_add(n + idx - 1, i);
                res = poly_add(res, poly_mul(tx_state->paks[i][j], _fec_inv(state, a_i)));
            }
        }
        out_pak[j] = res;
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

    for (i = 0; i < rx_state->num_info; i++) {
        y_i = rx_state->pak_xy_arr[i];

        fec_int_t pi_ycomp_y_div_ycomp_x_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, poly_add(y_i, rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, _fec_inv(state, poly_add(y_i, present_x[j])));
        }
        for (ii = 0; ii < pak_len; ii++) {
            rx_state->pak_arr[i][ii] = poly_mul(rx_state->pak_arr[i][ii], pi_ycomp_y_div_ycomp_x_i);
        }
    }

    for (i = 0; i < num_x_present; i++) {
        fec_int_t pi_xy_div_xx_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, poly_add(present_x[i], rx_state->missing_y[j]));
        }
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, _fec_inv(state, poly_add(present_x[i], present_x[j])));
        }

        unaligend_fec_int_t *pak = rx_state->pak_arr[n - num_x_present + i];
        for (ii = 0; ii < pak_len; ii++) {
            pak[ii] = poly_mul(pak[ii], pi_xy_div_xx_i);
        }
    }

    for (ii = 0; ii < pak_len; ii++) {
        fec_int_t res;

        for (i = 0; i < num_y_missing; i++) {

            res = 0;

            for (j = 0; j < rx_state->num_info; j++) {
                y_j = rx_state->pak_xy_arr[j];
                res = poly_add(res, poly_mul(rx_state->pak_arr[j][ii], _fec_inv(state, poly_add(y_j, rx_state->missing_y[i]))));
            }
            if (has_one_row) {
                res = poly_add(res, rx_state->ones_pak[ii]);
            }
            for(j = 0; j < num_x_present; j++) {
                res = poly_add(res, poly_mul(rx_state->pak_arr[n - num_x_present + j][ii], _fec_inv(state, poly_add(present_x[j], rx_state->missing_y[i]))));
            }

            rx_state->tmp_recovered_ints[i] = res;
        }

        for (i = 0; i < num_y_missing; i++) {
            if (i == 0 && has_one_row) {
                rx_state->ones_pak[ii] = rx_state->tmp_recovered_ints[i];
            } else {
                rx_state->pak_arr[n - num_x_present + i - has_one_row][ii] = rx_state->tmp_recovered_ints[i];
            }
        }
    }

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t pi_yx_div_yy_i = 1;
        for (j = 0; j < num_x_present; j++) {
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, poly_add(rx_state->missing_y[i], present_x[j]));
        }
        for (j = 0; j < num_y_missing; j++) {
            if(j == i) {
                continue;
            }
            pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, _fec_inv(state, poly_add(rx_state->missing_y[i], rx_state->missing_y[j])));
        }

        for (ii = 0; ii < pak_len; ii++) {
            if (i == 0 && has_one_row) {
                rx_state->ones_pak[ii] = poly_mul(rx_state->ones_pak[ii], pi_yx_div_yy_i);
            } else {
                rx_state->pak_arr[n - num_x_present + i - has_one_row][ii] = poly_mul(rx_state->pak_arr[n - num_x_present + i - has_one_row][ii], pi_yx_div_yy_i);
            }
        }
    }

    for (i = 0; i < rx_state->num_info; i++) {
        y_i = rx_state->pak_xy_arr[i];

        fec_int_t inv_pi_ycomp_y_div_ycomp_x_i = 1;
        for (j = 0; j < num_y_missing; j++) {
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i, _fec_inv(state, poly_add(y_i, rx_state->missing_y[j])));
        }
        for (j = 0; j < num_x_present; j++) {
            inv_pi_ycomp_y_div_ycomp_x_i = poly_mul(inv_pi_ycomp_y_div_ycomp_x_i,  poly_add(y_i, present_x[j]));
        }
        for (ii = 0; ii < pak_len; ii++) {
            rx_state->pak_arr[i][ii] = poly_mul(rx_state->pak_arr[i][ii], inv_pi_ycomp_y_div_ycomp_x_i);
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

#endif


void** fec_rx_get_info_paks(const fec_rx_state_t *rx_state) {
#ifndef FEC_LARGE_K
    return (void**)rx_state->info_paks;
#else
    return (void**)rx_state->pak_arr;
#endif
}

