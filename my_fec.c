
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define TRACE(fmt, ...) do { printf(fmt, ##__VA_ARGS__); } while(false)

#define VEC_TRACE(vec, num) \
{ \
    uint __i = 0;\
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
    fec_state_t* state;
    const unaligend_fec_int_t** paks;
} fec_tx_state_t;

typedef struct {
    fec_state_t* state;
    unaligend_fec_int_t** info_paks; // size = n
    unaligend_fec_int_t** redundancy_paks; // size = k
    fec_idx_t num_info;
    fec_idx_t num_redundant;

    fec_int_t *missing_y; // size = k
    fec_int_t *present_x; // size = k - 1

    fec_int_t *pi_xy_div_xx; // size = k - 1
    fec_int_t *pi_yx_div_yy; // size = k

    fec_int_t *present_y; // size = n - k

    fec_int_t *pi_ycomp_y_div_ycomp_x; // size = n - k

    // TODO: can replace with inplace in packet mem
    fec_int_t *tmp_vec_info; // size = n - k
    fec_int_t *tmp_vec_redundancy; // size = k - 1

    fec_int_t *tmp_recovered_ints; // size = k
    

} fec_rx_state_t;

// num must be > 0
uint log2_ceil(uint num) {
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

    // TODO: int overflows everywhere

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

    rx_state->state = state;
    CALLOC_ATTR(info_paks, state->n);

    CALLOC_ATTR(redundancy_paks, state->k);

    MALLOC_ATTR(tmp_vec_redundancy, state->k - 1);
    MALLOC_ATTR(missing_y, state->k);
    MALLOC_ATTR(present_x, state->k - 1);
    MALLOC_ATTR(pi_xy_div_xx, state->k - 1);
    MALLOC_ATTR(pi_yx_div_yy, state->k);
    MALLOC_ATTR(tmp_recovered_ints, state->k);
    MALLOC_ATTR(tmp_vec_info, state->n - state->k);
    MALLOC_ATTR(present_y, state->n - state->k);
    MALLOC_ATTR(pi_ycomp_y_div_ycomp_x, state->n - state->k);

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
    if (rx_state->info_paks != NULL) {
        free(rx_state->info_paks);
    }
    if (rx_state->redundancy_paks != NULL) {
        free(rx_state->redundancy_paks);
    }
    if (rx_state->tmp_vec_redundancy != NULL) {
        free(rx_state->tmp_vec_redundancy);
    }
    if (rx_state->missing_y != NULL) {
        free(rx_state->missing_y);
    }
    if (rx_state->present_x != NULL) {
        free(rx_state->present_x);
    }
    if (rx_state->pi_xy_div_xx != NULL) {
        free(rx_state->pi_xy_div_xx);
    }
    if (rx_state->pi_yx_div_yy != NULL) {
        free(rx_state->pi_yx_div_yy);
    }
    if (rx_state->tmp_recovered_ints != NULL) {
        free(rx_state->tmp_recovered_ints);
    }
    if (rx_state->tmp_vec_info != NULL) {
        free(rx_state->tmp_vec_info);
    }
    if (rx_state->present_y != NULL) {
        free(rx_state->present_y);
    }
    if (rx_state->pi_ycomp_y_div_ycomp_x != NULL) {
        free(rx_state->pi_ycomp_y_div_ycomp_x);
    }
}

bool fec_tx_add_info_pak(fec_tx_state_t *tx_state, const void* pak, fec_idx_t idx) {
    tx_state->paks[idx] = (const unaligend_fec_int_t*)pak;
    return true;
}

// TODO: idx here can overflow because n + k can be bigger than fec_int_t
bool fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_idx_t idx, bool *can_recover) {
    fec_int_t n = rx_state->state->n;
    if (idx < n) {
        rx_state->info_paks[idx] = (unaligend_fec_int_t*)pak;
        rx_state->num_info++;
    } else {
        rx_state->redundancy_paks[idx - n] = (unaligend_fec_int_t*)pak;
        rx_state->num_redundant++;
    }

    if (can_recover) {
        *can_recover = (rx_state->num_info + rx_state->num_redundant >= n);
    }
    
    return true;
}

// TODO: maybe
bool fec_tx_get_redundancy_pak(const fec_tx_state_t *tx_state, fec_idx_t idx, void *pak) {
    const fec_state_t *state = tx_state->state;
    fec_idx_t n = state->n;
    size_t pak_len = state->pak_len;
    fec_idx_t i;
    size_t j;
    unaligend_fec_int_t* out_pak = (unaligend_fec_int_t*)pak;

    for (i = 0; i < n; i++) {
        if (tx_state->paks[i] == NULL) {
            return false;
        }
    }

    // TODO: move if inside the fors?
    if (idx == 0) { 
        for (j = 0; j < pak_len; j++) {
            fec_int_t res = 0;
            for (i = 0; i < n; i++) {
                res = poly_add(res, tx_state->paks[i][j]);
            }
            out_pak[j] = res;
        }
    } else {
        for (j = 0; j < pak_len; j++) {
            fec_int_t res = 0;
            for (i = 0; i < n; i++) {
                fec_int_t a_i = poly_add(n + idx - 1, i);
                res = poly_add(res, poly_mul(tx_state->paks[i][j], _fec_inv(state, a_i)));
            }
            out_pak[j] = res;
        }
    }

    return true;
}

bool fec_rx_fill_missing_paks(fec_rx_state_t *rx_state) {
    const fec_state_t *state = rx_state->state;
    fec_idx_t n = state->n;
    fec_idx_t k = state->k;
    size_t pak_len = state->pak_len;

    fec_idx_t num_y_missing = n - rx_state->num_info;
    bool has_one_row;
    fec_idx_t num_x_present = 0;

    fec_idx_t i, j;
    fec_idx_t x_i, y_i, y_j;
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

    for (x_i = 0, i = 0; x_i < k - 1; x_i++) {
        if (rx_state->redundancy_paks[x_i + 1] == NULL) {
            continue;
        }
        rx_state->present_x[i] = n + x_i;
        i++;
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



int main(void) {

    fec_state_t state;
    fec_tx_state_t tx_state;
    fec_rx_state_t rx_state;

    uint i;

    uint16_t paks[3][1] = {{1}, {2}, {3}};
    uint16_t r_paks[3][1];

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

    for (i = 0; i < sizeof(paks)/sizeof(paks[0]); i++) {
        TRACE("%d ", rx_state.info_paks[i][0]);
    }
    TRACE("\n");

    fec_rx_reset(&rx_state); // not neede here, but to test

    fec_rx_destroy(&rx_state);
    fec_tx_destroy(&tx_state);
    fec_destroy(&state);

    return 0;
}
