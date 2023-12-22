
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

#define TRACE(fmt, ...) do { printf(fmt, ##__VA_ARGS__); } while(false)

#define VEC_TRACE(vec, num) \
{ \
    unsigned int __i = 0;\
    TRACE("%s:", #vec);\
    for(__i = 0; __i < (num); __i++) {\
        TRACE(" %d", (vec)[__i]);\
    }\
    TRACE("\n");\
}

#define CHECK(cond) if(!(cond)) { TRACE("check failed %d %s\n", __LINE__, #cond); }

typedef uint16_t fec_int_t;
typedef uint32_t fec_idx_t;

typedef fec_int_t __attribute__((aligned(1))) unaligend_fec_int_t;

typedef struct {
    fec_idx_t n;
    fec_idx_t k;
    size_t pak_len;

    // n+k-1 <= (1<<(sizeof(fec_int_t)*8))

    // 0,...,n+k-2 are the numbers
    // 1,...,1<<ceil(log2(n+k-1))

    fec_int_t* inv_arr; // upper power of two closest to n+k-2

} fec_state_t;

#define POLY_G ((fec_int_t)0b0000000000101011)

typedef struct {
    const fec_state_t* state;
    const unaligend_fec_int_t** paks;
} fec_tx_state_t;

typedef struct {
    const fec_state_t* state;
    unaligend_fec_int_t** info_paks; // size = n
    unaligend_fec_int_t** redundancy_paks; // size = real k
    fec_idx_t num_info;
    fec_idx_t num_redundant;

    // all use this:
    fec_int_t *missing_y; // size = k
    fec_int_t *present_x; // size = k - 1

#ifdef FEC_MIN_MEM
    // min mem uses this:
    fec_int_t *tmp_recovered_ints; // size = k
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

// num must be > 0
unsigned int log2_ceil(unsigned int num) {
    if (num == 1) {
        return 0;
    }
    return (sizeof(num)*8) - __builtin_clz(num - 1);
}

fec_int_t poly_add(fec_int_t a, fec_int_t b) {
    return a ^ b;
}

fec_int_t poly_shift_1(fec_int_t a) {
    if ((a & (1<<(sizeof(fec_int_t)*8 - 1))) != 0) {
        return poly_add((a << 1), POLY_G);
    } else {
        return (a << 1);
    }
}

fec_int_t poly_mul(fec_int_t a, fec_int_t b) {
    fec_int_t cur_bit = (1 << (sizeof(fec_int_t)*8 - 1));
    size_t i;
    fec_int_t res = 0;

    for (i = 0; i < sizeof(fec_int_t)*8; i++) {
        res = poly_shift_1(res);
        if ((b & cur_bit) != 0) {
            res = poly_add(res, a);
        }
        cur_bit >>= 1;
    }

    return res;
}

fec_int_t poly_pow(fec_int_t a, fec_int_t n) {
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

fec_int_t poly_inv(fec_int_t a) {
    fec_int_t res = poly_pow(a, (fec_int_t)-2);
    CHECK(poly_mul(a, res) == 1);
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

void fec_rx_destroy(fec_rx_state_t *rx_state);

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
    CALLOC_ATTR(info_paks, n);
    CALLOC_ATTR(redundancy_paks, k);
    MALLOC_ATTR(missing_y, MIN(k, n));
    MALLOC_ATTR(present_x, MIN(k - 1, n));

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
    return true;
}

void fec_rx_reset(fec_rx_state_t *rx_state) {
    memset(rx_state->info_paks, 0, rx_state->state->n * sizeof(rx_state->info_paks[0]));
    memset(rx_state->redundancy_paks, 0, rx_state->state->k * sizeof(rx_state->redundancy_paks[0]));
    rx_state->num_info = 0;
    rx_state->num_redundant = 0;
}

void fec_rx_destroy(fec_rx_state_t *rx_state) {
#define FREE_ATTR(name) \
    if (rx_state->name != NULL) { \
        free(rx_state->name); \
    }
    
    FREE_ATTR(info_paks);
    FREE_ATTR(redundancy_paks);
    FREE_ATTR(present_x);
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

bool fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_idx_t idx, bool *can_recover) {
    fec_int_t n = rx_state->state->n;

    if (idx >= n + rx_state->state->k) {
        return false;
    }

    if (idx < n) {
        if (rx_state->info_paks[idx] == NULL) {
            rx_state->info_paks[idx] = (unaligend_fec_int_t*)pak;
            rx_state->num_info++;
        }
    } else {
        if (rx_state->redundancy_paks[idx - n] == NULL) {
            rx_state->redundancy_paks[idx - n] = (unaligend_fec_int_t*)pak;
            if (idx != n) {
                bool has_one_row = (rx_state->redundancy_paks[0] != NULL);
                rx_state->present_x[rx_state->num_redundant - has_one_row] = idx - 1;
            }
            rx_state->num_redundant++;
        }
    }

    if (can_recover) {
        *can_recover = (rx_state->num_info + rx_state->num_redundant >= n);
    }
    
    return true;
}

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
        fec_int_t tmp_vec_1s;
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
#endif

int main(void) {

    fec_state_t state;
    fec_tx_state_t tx_state;
    fec_rx_state_t rx_state;

    unsigned int i, j;

    uint16_t paks[3][2] = {{1, 4}, {2, 5}, {3, 6}};
    uint16_t r_paks[3][2];

    TRACE("--0--\n");

    fec_init(&state, sizeof(paks)/sizeof(paks[0]), sizeof(r_paks)/sizeof(r_paks[0]), sizeof(paks[0])/sizeof(paks[0][0]));
    fec_tx_init(&tx_state, &state);
    fec_rx_init(&rx_state, &state);

    TRACE("--1--\n");

    for (i = 0; i < sizeof(paks)/sizeof(paks[0]); i++) {
        fec_tx_add_info_pak(&tx_state, paks[i], i);
    }

    TRACE("--2--\n");

    for (i = 0; i < sizeof(r_paks)/sizeof(r_paks[0]); i++) {
        fec_tx_get_redundancy_pak(&tx_state, i, r_paks[i]);
    }

    TRACE("--3--\n");

    for (i = 0; i < sizeof(paks)/sizeof(paks[0]); i++) {
        if (i == 0 || i == 2 || i == 1) {
            continue;
        }
        fec_rx_add_pak(&rx_state, paks[i], i, NULL);
    }

    TRACE("--4--\n");

    for (i = 0; i < sizeof(r_paks)/sizeof(r_paks[0]); i++) {
        fec_rx_add_pak(&rx_state, r_paks[i], (sizeof(paks)/sizeof(paks[0])) + i, NULL);
    }

    TRACE("--5--\n");

    if (!fec_rx_fill_missing_paks(&rx_state)) {
        return 0;
    }

    TRACE("--6--\n");

    for (j = 0; j < sizeof(paks[0])/sizeof(paks[0][0]); j++) {
        for (i = 0; i < sizeof(paks)/sizeof(paks[0]); i++) {
            TRACE("%d ", rx_state.info_paks[i][j]);
        }
        TRACE("\n");
    }

    fec_rx_reset(&rx_state); // not neede here, but to test

    fec_rx_destroy(&rx_state);
    fec_tx_destroy(&tx_state);
    fec_destroy(&state);

    return 0;
}
