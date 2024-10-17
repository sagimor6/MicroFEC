#ifndef WIN32
#define _GNU_SOURCE
#endif
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

uint32_t my_hash(const void* buf, size_t size) {
    size_t i;
    uint32_t res = 0;
    for (i = 0; i < size; i++) {
        res += ((const uint8_t*)buf)[i] * (i + 1);
    }

    return res;
}

void rand_fill(void* buf, size_t len) {
    size_t i;
    for (i = 0; i < len; i++) {
        // TODO: inefficient but idc
        ((uint8_t*)buf)[i] = (uint8_t)rand_range(1<<8);
    }
}

#ifdef _WIN32
#include <realtimeapiset.h>
#include <processthreadsapi.h>
#include <profileapi.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
uint64_t get_timestamp() {
    // uint64_t t;
    // CHECK(QueryThreadCycleTime(GetCurrentThread(), &t));
    // t /= 4;

    uint64_t t, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&t);
    return (uint64_t)(((double)t * 1000000000) / freq);
cleanup:
    return t;
}
#else
#include <sched.h>
uint64_t get_timestamp() {
    struct timespec tp = {0};
    //CHECK(clock_gettime(CLOCK_MONOTONIC, &tp) == 0);
    CHECK(clock_gettime(CLOCK_MONOTONIC_RAW, &tp) == 0);

cleanup:
    return (tp.tv_sec*1000000000ULL) + tp.tv_nsec;
}
#endif

bool test_perf(unsigned int n, unsigned int k, unsigned int pak_len) {
    // const unsigned int n = 1000;
    // const unsigned int k = 100;
    // const size_t pak_len = 1400 / sizeof(uint16_t);
    fec_int_t *paks = NULL;
    fec_int_t *r_paks = NULL;
    unsigned int *rcv_idxs = NULL;
    unsigned int i;
    bool inited_inv_cache = false;
    fec_inv_cache_t inv_cache;
    bool inited_tx_state = false;
    fec_tx_state_t tx_state;
    bool inited_rx_state = false;
    fec_rx_state_t rx_state;
    uint64_t start_time, end_time;
    fec_int_t *rx_dest_buf = NULL;
    bool success = false;

#ifdef _WIN32
    // we want accuracy
    CHECK(SetThreadAffinityMask(GetCurrentThread(), 1) != 0);
    CHECK(SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS));
    CHECK(SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL));
#else
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(0, &cpu_set);
    
    CHECK(sched_setaffinity(0, sizeof(cpu_set), &cpu_set) == 0);

    struct sched_param sched_param = {0};
    sched_param.sched_priority = sched_get_priority_max(SCHED_FIFO);
    CHECK(sched_param.sched_priority >= 0);
    CHECK(sched_setscheduler(0, SCHED_FIFO, &sched_param) == 0);
#endif

    CHECK((pak_len % sizeof(fec_int_t)) == 0);

    pak_len /= sizeof(fec_int_t);

    paks = malloc(n * pak_len * sizeof(uint16_t));
    CHECK(paks != NULL);

    r_paks = malloc(k * pak_len * sizeof(uint16_t));;
    CHECK(r_paks != NULL);

    rx_dest_buf = malloc(n * pak_len * sizeof(uint16_t));
    CHECK(rx_dest_buf != NULL);

    srand(1337/*time(NULL)*/);

    rcv_idxs = rand_unique_indexes(n + k, n);
    CHECK(rcv_idxs != NULL);

    //rand_fill(paks, n * pak_len * sizeof(uint16_t));

    {
        for(i = 0; i < n * pak_len; i++) {
            paks[i] = i;
        }
    }

#ifndef FEC_USER_GIVEN_BUFFER
    memcpy(rx_dest_buf, paks, n * pak_len * sizeof(uint16_t));
#endif

#ifdef FEC_HAS_CLMUL64
    TRACE("FEC_HAS_CLMUL64 ");
#endif
#ifdef FEC_HAS_CLMUL32
    TRACE("FEC_HAS_CLMUL32 ");
#endif
#ifdef FEC_HAS_128_INT_VEC
    TRACE("FEC_HAS_128_INT_VEC ");
#endif
#ifdef FEC_HAS_64_INT_VEC
    TRACE("FEC_HAS_64_INT_VEC ");
#endif
#ifdef FEC_HAS_64BIT
    TRACE("FEC_HAS_64BIT ");
#endif
#ifdef FEC_HAS_32BIT
    TRACE("FEC_HAS_32BIT ");
#endif
    TRACE("\n");

    TRACE("--0--\n");
    start_time = get_timestamp();

    CHECK(fec_inv_cache_init(&inv_cache, n, k) == FEC_STATUS_SUCCESS);
    inited_inv_cache = true;
    CHECK(fec_tx_init(&tx_state, n, pak_len) == FEC_STATUS_SUCCESS);
    inited_tx_state = true;
#ifndef FEC_USER_GIVEN_BUFFER
    CHECK(fec_rx_init(&rx_state, n, k, pak_len) == FEC_STATUS_SUCCESS);
#else
    CHECK(fec_rx_init(&rx_state, n, k, pak_len, rx_dest_buf) == FEC_STATUS_SUCCESS);
#endif
    inited_rx_state = true;

    PRINT_TS_DIFF("--total init time--");

    for (i = 0; i < n; i++) {
        CHECK(fec_tx_add_info_pak(&tx_state, &paks[i*pak_len], i) ==  FEC_STATUS_SUCCESS);
    }

    PRINT_TS_DIFF("--fec_tx_add_info_pak--");

    for (i = 0; i < k; i++) {
        CHECK(fec_tx_get_redundancy_pak(&tx_state, &inv_cache, i, &r_paks[i*pak_len]) == FEC_STATUS_SUCCESS);
    }

    PRINT_TS_DIFF("--fec_tx_get_redundancy_pak--");

    for (i = 0; i < n; i++) {
        unsigned int idx = rcv_idxs[i];
        fec_status_t status;
        
        if (idx < n) {
            status = fec_rx_add_pak(&rx_state, &paks[idx*pak_len], idx);
        } else {
            status = fec_rx_add_pak(&rx_state, &r_paks[(idx - n)*pak_len], idx);
        }
        if (i == n - 1) {
            CHECK(status == FEC_STATUS_SUCCESS);
        } else {
            CHECK(status == FEC_STATUS_MORE_PACKETS_NEEDED);
        }
        
        if (idx < n) {
            status = fec_rx_add_pak(&rx_state, &paks[idx*pak_len], idx);
        } else {
            status = fec_rx_add_pak(&rx_state, &r_paks[(idx - n)*pak_len], idx);
        }
        if (i == n - 1) {
            CHECK(status == FEC_STATUS_CAN_DROP_ALREADY_RECOVERABLE);
        } else {
            CHECK(status == FEC_STATUS_CAN_DROP_DUP_PAK);
        }
    }

    PRINT_TS_DIFF("--fec_rx_add_pak--");

    CHECK(fec_rx_fill_missing_paks(&rx_state, &inv_cache) == FEC_STATUS_SUCCESS);

    PRINT_TS_DIFF("--fec_rx_fill_missing_paks--");

#ifndef FEC_USER_GIVEN_BUFFER
    uint16_t** res = (uint16_t**)fec_rx_get_info_paks(&rx_state);
    for (i = 0; i < n; i++) {
        CHECK(memcmp(res[i], &rx_dest_buf[i*pak_len], pak_len * sizeof(uint16_t)) == 0);
    }
#else
    CHECK(memcmp(rx_dest_buf, paks, n * pak_len * sizeof(uint16_t)) == 0);
#endif

    success = true;

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
    if (rx_dest_buf != NULL) {
        free(rx_dest_buf);
    }

    if (inited_rx_state) {
        fec_rx_destroy(&rx_state);
    }

    if (inited_tx_state) {
        fec_tx_destroy(&tx_state);
    }

    if (inited_inv_cache) {
        fec_inv_cache_destroy(&inv_cache);
    }

    return success;
}

int main(int argc, const char* argv[]) {
    if (argc != 4) {
        TRACE("usage: %s N K PAK_LEN\n", argc > 0 ? argv[0] : "./program");
        return 0;
    }

    bool success = test_perf(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    return success ? 0 : 1;
}
