#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "micro_fec.h"

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

#define CHECK(cond) if(!(cond)) { TRACE("check failed %d %s\n", __LINE__, #cond); goto cleanup; }


unsigned int rand_range(unsigned int n) {
    if (n > ((unsigned int)RAND_MAX + 1)) {
        abort();
    }

    unsigned int max = RAND_MAX - (((unsigned int)RAND_MAX + 1) % n);

    while (true) {
        unsigned int val = rand();
        if (val <= max) {
            return (val % n);
        }
    }
}

// TODO: this is can be done with O(k) mem with dict.
unsigned int* rand_unique_indexes(unsigned int n, unsigned int k) {
    unsigned int* arr = malloc(n*sizeof(unsigned int));
    unsigned int i;

    if (arr == NULL) {
        return NULL;
    }

    for (i = 0; i < n; i++) {
        arr[i] = i;
    }

    for (i = 0; i < k; i++) {
        unsigned int j = rand_range(n - i) + i;
        unsigned int tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }

    unsigned int* arr2 = realloc(arr, k*sizeof(unsigned int));
    if (arr2 == NULL) {
        return arr; // weird but ok
    }

    return arr2;
}

void rand_fill(void* buf, size_t len) {
    size_t i;
    for (i = 0; i < len; i++) {
        // TODO: inefficient but idc
        ((uint8_t*)buf)[i] = (uint8_t)rand_range(1<<8);
    }
}

uint64_t get_timestamp() {
    struct timespec tp = {0};
    //CHECK(clock_gettime(CLOCK_MONOTONIC, &tp) == 0);
    CHECK(clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp) == 0);

cleanup:
    return (tp.tv_sec*1000000000ULL) + tp.tv_nsec;
}

void test_perf(unsigned int n, unsigned int k, unsigned int pak_len) {
    // const unsigned int n = 1000;
    // const unsigned int k = 100;
    // const size_t pak_len = 1400 / sizeof(uint16_t);
    fec_int_t *paks = NULL;
    fec_int_t *r_paks = NULL;
    unsigned int *rcv_idxs = NULL;
    unsigned int i;
    bool inited_state = false;
    fec_state_t state;
    bool inited_tx_state = false;
    fec_tx_state_t tx_state;
    bool inited_rx_state = false;
    fec_rx_state_t rx_state;
    uint64_t start_time, end_time;

    CHECK((pak_len % sizeof(fec_int_t)) == 0);

    pak_len /= sizeof(fec_int_t);

    paks = malloc(n * pak_len * sizeof(uint16_t));
    CHECK(paks != NULL);

    r_paks = malloc(k * pak_len * sizeof(uint16_t));;
    CHECK(r_paks != NULL);

    srand(1337/*time(NULL)*/);

    rcv_idxs = rand_unique_indexes(n + k, n);
    CHECK(rcv_idxs != NULL);

    rand_fill(paks, n * pak_len * sizeof(uint16_t));

    TRACE("--0--\n");
    start_time = get_timestamp();

    CHECK(fec_init(&state, n, k, pak_len));
    inited_state = true;
    CHECK(fec_tx_init(&tx_state, &state));
    inited_tx_state = true;
    CHECK(fec_rx_init(&rx_state, &state));
    inited_rx_state = true;

    end_time = get_timestamp();
    TRACE("--total init time-- %f\n", (end_time - start_time) / ((double)1000000000));
    start_time = get_timestamp();

    for (i = 0; i < n; i++) {
        fec_tx_add_info_pak(&tx_state, &paks[i*pak_len], i);
    }

    end_time = get_timestamp();
    TRACE("--fec_tx_add_info_pak-- %f\n", (end_time - start_time) / ((double)1000000000));
    start_time = get_timestamp();

    for (i = 0; i < k; i++) {
        fec_tx_get_redundancy_pak(&tx_state, i, &r_paks[i*pak_len]);
    }

    end_time = get_timestamp();
    TRACE("--fec_tx_get_redundancy_pak-- %f\n", (end_time - start_time) / ((double)1000000000));
    start_time = get_timestamp();

    for (i = 0; i < n; i++) {
        unsigned int idx = rcv_idxs[i];
        bool can_recover;
        bool can_discard;
        
        if (idx < n) {
            CHECK(fec_rx_add_pak(&rx_state, &paks[idx*pak_len], idx, &can_recover, &can_discard));
        } else {
            CHECK(fec_rx_add_pak(&rx_state, &r_paks[(idx - n)*pak_len], idx, &can_recover, &can_discard));
        }
        CHECK(can_recover == (i == n - 1));
        CHECK(!can_discard);
        
        if (idx < n) {
            CHECK(fec_rx_add_pak(&rx_state, &paks[idx*pak_len], idx, &can_recover, &can_discard));
        } else {
            CHECK(fec_rx_add_pak(&rx_state, &r_paks[(idx - n)*pak_len], idx, &can_recover, &can_discard));
        }
        CHECK(can_recover == (i == n - 1));
        CHECK(can_discard);
    }

    end_time = get_timestamp();
    TRACE("--fec_rx_add_pak-- %f\n", (end_time - start_time) / ((double)1000000000));
    start_time = get_timestamp();

    CHECK(fec_rx_fill_missing_paks(&rx_state));

    end_time = get_timestamp();
    TRACE("--fec_rx_fill_missing_paks-- %f\n", (end_time - start_time) / ((double)1000000000));
    start_time = get_timestamp();

    uint16_t** res = (uint16_t**)fec_rx_get_info_paks(&rx_state);

    for (i = 0; i < n; i++) {
        if (res[i] == &paks[i*pak_len]) {
            continue;
        }
        CHECK(memcmp(res[i], &paks[i*pak_len], pak_len * sizeof(uint16_t)) == 0);
    }

cleanup:
    if (paks != NULL) {
        free(paks);
    }
    if (r_paks != NULL) {
        free(r_paks);
    }
    if (rcv_idxs != NULL) {
        free(rcv_idxs);
    }

    if (inited_rx_state) {
        fec_rx_destroy(&rx_state);
    }

    if (inited_tx_state) {
        fec_tx_destroy(&tx_state);
    }

    if (inited_state) {
        fec_destroy(&state);
    }
}

int main(int argc, const char* argv[]) {
    if (argc != 4) {
        TRACE("usage: %s N K PAK_LEN\n", argc > 0 ? argv[0] : "./program");
        return 0;
    }

    test_perf(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    return 0;
}
