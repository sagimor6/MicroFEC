
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>

typedef uint16_t fec_int_t;

typedef fec_int_t __attribute__((aligned(1))) unaligend_fec_int_t;

typedef struct {
    fec_int_t n;
    fec_int_t k;
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
    fec_int_t num_info;
    fec_int_t num_redundant;

    fec_int_t *missing_y; // size = k
    fec_int_t *present_x; // size = k
    bool missing_1s;

    fec_int_t *pi_xx; // size = k
    fec_int_t *pi_yy; // size = k
    fec_int_t *pi_xy; // size = k
    fec_int_t *pi_yx; // size = k

    fec_int_t *present_y; // size = n - k
    fec_int_t *pi_ycomp_y; // size = n - k
    fec_int_t *pi_ycomp_x; // size = n - k

    

    // TODO: can replace with inplace in packet mem
    fec_int_t *tmp_vec_info; // size = n - k
    fec_int_t tmp_vec_1s;
    fec_int_t *tmp_vec_redundancy; // size = k

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
    return poly_pow(a, (fec_int_t)-2);
}

bool fec_init(fec_state_t *state, fec_int_t n, fec_int_t k, size_t pak_len) {
    state->n = n;
    state->k = k;

    // TODO: int overflows everywhere

    if (n + k - 1 > 0) {
        size_t arr_size = (1 << log2_ceil(n + k - 1)) - 1;
        state->inv_arr = malloc(sizeof(fec_int_t) * arr_size);
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

bool fec_tx_init(fec_tx_state_t *tx_state, fec_state_t *state) {
    tx_state->state = state;
    tx_state->paks = calloc(state->n, sizeof(tx_state->paks[0]));
    if (tx_state->paks == NULL) {
        return false;
    }

    return true;
}

bool fec_rx_init(fec_rx_state_t *rx_state, fec_state_t *state) {
    rx_state->state = state;
    rx_state->info_paks = calloc(state->n, sizeof(rx_state->info_paks[0]));
    if (rx_state->info_paks == NULL) {
        return false;
    }

    return true;
}

bool fec_tx_add_info_pak(fec_tx_state_t *tx_state, const void* pak, fec_int_t idx) {
    tx_state->paks[idx] = (const unaligend_fec_int_t*)pak;
    return true;
}

// TODO: idx here can overflow because n + k can be bigger than fec_int_t
bool fec_rx_add_pak(fec_rx_state_t *rx_state, void* pak, fec_int_t idx, bool *can_recover) {
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
bool fec_tx_get_redundancy_pak(const fec_tx_state_t *tx_state, fec_int_t idx, void *pak) {
    const fec_state_t *state = tx_state->state;
    fec_int_t n = state->n;
    fec_int_t k = state->k;
    size_t pak_len = state->pak_len;
    fec_int_t i;
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
                res = poly_add(res, poly_mul(tx_state->paks[i][j], state->inv_arr[a_i - 1]));
            }
            out_pak[j] = res;
        }
    }

    return true;
}

bool fec_rx_fill_missing_paks(fec_rx_state_t *rx_state) {
    const fec_state_t *state = rx_state->state;
    fec_int_t n = state->n;
    fec_int_t k = state->k;
    size_t pak_len = state->pak_len;

    fec_int_t num_y_missing = n - rx_state->num_info;
    fec_int_t num_x_present = rx_state->num_redundant;

    fec_int_t i, j;
    fec_int_t idx_i, idx_j;
    fec_int_t x_i, y_i, y_j, x_j;
    size_t ii;

    for (y_i = 0, i = 0, j = 0; y_i < n; y_i++) {
        if (rx_state->info_paks[y_i] != NULL) {
            rx_state->present_y[j] = y_i;
            j++;
        } else {
            rx_state->missing_y[i] = y_i;
            i++;
        }
    }

    rx_state->missing_1s = (rx_state->redundancy_paks[0] == NULL);

    for (x_i = 0, i = 0; x_i < k - 1; x_i++) {
        if (rx_state->redundancy_paks[x_i + 1] == NULL) {
            continue;
        }
        rx_state->present_x[i] = n + x_i;
        i++;
    }



    for (i = 0, idx_i = 0, idx_j = 0; i < k - 1; i++) {
        if (rx_state->redundancy_paks[i+1] != NULL) {
            continue;
        }

        fec_int_t res = 1;

        for (j = 0; j < k - 1; j++) {
            if (rx_state->redundancy_paks[j+1] != NULL) {
                continue;
            }
            if(i == j) {
                idx_j++;
                continue;
            }

            res = poly_mul(res, state->inv_arr[poly_add(n + i, n + j)]);

            idx_j++;
        }

        rx_state->pi_xx[i] = res; // TODO: i or idx_i?

        idx_i++;
    }

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_y_missing; j++) {
            if(j == i) {
                continue;
            }
            res = poly_mul(res, state->inv_arr[poly_add(rx_state->missing_y[i], rx_state->missing_y[j])]);
        }
        rx_state->pi_yy[i] = res;
    }

    for (i = 0; i < num_x_present; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_x_present; j++) {
            if(j == i) {
                continue;
            }
            res = poly_mul(res, state->inv_arr[poly_add(rx_state->present_x[i], rx_state->present_x[j])]);
        }
        rx_state->pi_xx[i] = res;
    }

    for (i = 0; i < num_x_present; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_y_missing; j++) {
            res = poly_mul(res, poly_add(rx_state->present_x[i], rx_state->missing_y[j]));
        }
        rx_state->pi_xy[i] = res;
    }

    for (i = 0; i < num_y_missing; i++) {
        fec_int_t res = 1;
        for (j = 0; j < num_x_present; j++) {
            res = poly_mul(res, poly_add(rx_state->missing_y[i], rx_state->present_x[j]));
        }
        rx_state->pi_yx[i] = res;
    }

    for (y_i = 0, i = 0; y_i < n; y_i++) {
        if (rx_state->info_paks[y_i] == NULL) {
            continue;
        }
        fec_int_t res = 1;
        for (j = 0; j < num_y_missing; j++) {
            res = poly_mul(res, poly_add(y_i, rx_state->missing_y[j]));
        }
        rx_state->pi_ycomp_y[i] = res;
        i++;
    }

    for (y_i = 0, i = 0; y_i < n; y_i++) {
        if (rx_state->info_paks[y_i] == NULL) {
            continue;
        }
        fec_int_t res = 1;
        for (j = 0; j < num_x_present; j++) {
            res = poly_mul(res, state->inv_arr[poly_add(y_i, rx_state->present_x[j])]);
        }
        rx_state->pi_ycomp_x[i] = res;
        i++;
    }

    for (ii = 0; ii < pak_len; ii++) {
        fec_int_t res;
        fec_int_t tmp_vec_1s;
        for (y_i = 0, i = 0; y_i < n; y_i++) {
            if (rx_state->info_paks[y_i] == NULL) {
                continue;
            }
            res = poly_mul(rx_state->info_paks[y_i][ii], rx_state->pi_ycomp_y[i]);
            res = poly_mul(res, rx_state->pi_ycomp_x[i]);
            rx_state->tmp_vec_info[i] = res;
            i++;
        }
        if(rx_state->redundancy_paks[0] != NULL) {
            tmp_vec_1s = rx_state->redundancy_paks[0][ii];
        }
        for (x_i = 0, i = 0; x_i < k - 1; x_i++) {
            if (rx_state->redundancy_paks[x_i] == NULL) {
                continue;
            }
            res = poly_mul(rx_state->redundancy_paks[x_i + 1][ii], rx_state->pi_xy[i]);
            res = poly_mul(res, rx_state->pi_xx[i]);
            rx_state->tmp_vec_redundancy[i] = res;
            i++;
        }

        for (i = 0; i < k; i++) {

            res = 0;

            for (y_j = 0, j = 0; y_j < n; y_j++) {
                if (rx_state->info_paks[y_j] == NULL) {
                    continue;
                }
                res = poly_add(res, poly_mul(rx_state->tmp_vec_info[j], state->inv_arr[poly_add(y_j, rx_state->missing_y[i])]));
                j++;
            }
            if (rx_state->redundancy_paks[0] != NULL) {
                res = poly_add(res, tmp_vec_1s);
            }
            for (x_j = 0, j = 0; x_j < k - 1; x_j++) {
                if (rx_state->redundancy_paks[x_j] == NULL) {
                    continue;
                }
                res = poly_add(res, poly_mul(rx_state->tmp_vec_redundancy[j], state->inv_arr[poly_add(x_j, rx_state->missing_y[i])]));
                j++;
            }

            res = poly_mul(res, rx_state->pi_yx[i]);
            res = poly_mul(res, rx_state->pi_yy[i]);

            rx_state->redundancy_paks[rx_state->present_x[i] - n][ii] = res;
        }


    }

}



int main(void) {
    return 0;
}
