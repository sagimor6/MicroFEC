#include <stdio.h>

#include "my_fec.h"

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

    uint16_t** res = (uint16_t**)fec_rx_get_info_paks(&rx_state);

    for (j = 0; j < sizeof(paks[0])/sizeof(paks[0][0]); j++) {
        for (i = 0; i < sizeof(paks)/sizeof(paks[0]); i++) {
            TRACE("%d ", res[i][j]);
        }
        TRACE("\n");
    }

    fec_rx_reset(&rx_state); // not neede here, but to test

    fec_rx_destroy(&rx_state);
    fec_tx_destroy(&tx_state);
    fec_destroy(&state);

    return 0;
}
