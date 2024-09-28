#include <stdint.h>

#include "clmul.h"

#include "dummy_defs.h"

#if INIT_FUNC_NAME != -1
void INIT_FUNC_NAME(fec_perf_int_t* restrict perf_col, fec_int_t val, LEN_PARAM_TYPE len) {
    LEN_PARAM_TYPE i;
    union {
        fec_int_t val_arr[sizeof(fec_perf_int_t)/sizeof(fec_int_t)];
        fec_perf_int_t vec_val;
    } vec_val = {.val_arr = {
#if defined (FEC_HAS_CLMUL32) || defined(FEC_HAS_64_INT_VEC) || !defined(FEC_HAS_64BIT)
        [0] = val
#else
        [15] = val
#endif
    }};
    for(i = 0; i < len; i++) {
        ((__typeof__(vec_val)*)perf_col)[i] = vec_val;
    }
}
#endif

void PERF_DEBUG_ATTRS FMA_FUNC_NAME(fec_perf_int_t* restrict perf_col, LEN_PARAM_TYPE len, fec_int_t a, INPUT_ARGS) {
    LEN_PARAM_TYPE i;

#if defined(FEC_HAS_CLMUL64)
    _poly_t _a = _POLY_1VAL(a);

    for (i = 0; i < len - 1; i += 2) {
        _poly_t _b = _POLY_2VAL(READ_INPUT(i), READ_INPUT(i+1));

        _poly_t _c = _POLY_CLMUL(_a, _b);

        perf_col[i] ^= _POLY_EXTRACT(_c, uint32_t, 0);
        perf_col[i+1] ^= _POLY_EXTRACT(_c, uint32_t, 1);
    }
    
    if (i == len - 1) {
        _poly_t _b = _POLY_1VAL(READ_INPUT(i));

        _poly_t _c = _POLY_CLMUL(_a, _b);

        perf_col[i] ^= _POLY_EXTRACT(_c, uint32_t, 0);
    }
#elif defined(FEC_HAS_CLMUL32)
    _poly_t _a = _POLY_1VAL(a);

    for (i = 0; i < len; i++) {
        _poly_t _b = _POLY_1VAL(READ_INPUT(i));

        _poly_t _c = _POLY_CLMUL(_a, _b);

        perf_col[i] ^= _POLY_EXTRACT(_c, uint32_t, 0);
    }
#elif defined(FEC_HAS_128_INT_VEC) || (defined(FEC_HAS_64_INT_VEC) && !(defined(__x86_64__) || defined(__i386__)))
    u16x16 shifts = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    u16x16 _a = (u16x16)(((i16x16)(a << shifts)) >> 15);

    for (i = 0; i < len; i++) {
        u16x16 _c = _a & READ_INPUT(i);
        perf_col[i] ^= _c;
    }
#elif defined(FEC_HAS_64_INT_VEC) && (defined(__x86_64__) || defined(__i386__)) && defined(__SSE__)
    int j;
    union {
        __m64 mmx[4];
        __m128 sse[2];
    } _a_tmp;
    const __m64 __a = _mm_set1_pi16(a);
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        __m64 res = _mm_srai_pi16(_mm_mullo_pi16(__a, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}), 15);
        asm inline(
                "movq %1, %0;"
            : "=m"(_a_tmp.mmx[j])
            : "y"(res)
            : );
    }

     __m128 _a_low;
     __m128 _a_high;

    asm inline(
                "movaps %1, %0;"
            : "=x"(_a_low)
            : "m"(_a_tmp.sse[0])
            : );
    
    asm inline(
                "movaps %1, %0;"
            : "=x"(_a_high)
            : "m"(_a_tmp.sse[1])
            : );

    for (i = 0; i < len; i++) {

        uint16_t pak_val = READ_INPUT(i);

        uint32_t pak_val_2 =  (pak_val | (((uint32_t)pak_val) << 16));
        __m128 pak_val_8;

        __m128 res;

        asm inline(
                "movss %[pak_val_2], %[pak_val_8];"
                "shufps $0, %[pak_val_8], %[pak_val_8];"
            : [pak_val_8] "=x"(pak_val_8)
            : [pak_val_2] "m"(pak_val_2)
            : );

        res = pak_val_8;

        asm inline(
                "andps %[_a_low], %[res];"
            : [res] "+x"(res)
            : [_a_low] "x"(_a_low)
            : );

        asm inline(
                "xorps %0, %1;"
                "movaps %1, %0;"
            : "+m"(perf_col[i][0]), "+x"(res)
            : 
            : );

        res = pak_val_8;

        asm inline(
                "andps %[_a_high], %[res];"
            : [res] "+x"(res)
            : [_a_high] "x"(_a_high)
            : );

        asm inline(
                "xorps %0, %1;"
                "movaps %1, %0;"
            : "+m"(perf_col[i][1]), "+x"(res)
            : 
            : );

        // perf_col[i][0] = _mm_xor_ps(perf_col[i][0], _mm_and_ps(_a_low, pak_val_8));
        // perf_col[i][1] = _mm_xor_ps(perf_col[i][1], _mm_and_ps(_a_high, pak_val_8));
    }
#elif defined(FEC_HAS_64_INT_VEC) && (defined(__x86_64__) || defined(__i386__)) && defined(__MMX__)
    int j;
    __m64 _a[4];
    const __m64 __a = _mm_set1_pi16(a);
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        _a[j] = _mm_srai_pi16(_mm_mullo_pi16(__a, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}), 15);
        //_a[j] = _mm_cmpgt_pi16((__m64){0}, _mm_mullo_pi16(__a, (__m64)(u16x4){1 << (15 - (j*4 + 0)), 1 << (15 - (j*4 + 1)), 1 << (15 - (j*4 + 2)), 1 << (15 - (j*4 + 3))}));
    }

    for (i = 0; i < len; i++) {
        __m64 pak_val = _mm_set1_pi16(READ_INPUT(i));
        #pragma GCC unroll 4
        for (j = 0; j < 4; j++) {

            // perf_col[i][j] = _mm_xor_si64(perf_col[i][j], _mm_and_si64(_a[j], pak_val));

            __m64 res = _mm_and_si64(_a[j], pak_val);

            asm inline(
                "pxor %0, %1;"
                "movq %1, %0;"
            : "+m"(perf_col[i][j]), "+y"(res)
            : 
            : );
        }
    }
#elif defined(FEC_HAS_64BIT)
    int j;
    uint64_t _a[4];
    
    #pragma GCC unroll 4
    for (j = 0; j < 4; j++) {
        _a[j] = a * ((1ULL<<((16+1)*(4*j+0) - 64*j)) + (1ULL<<((16+1)*(4*j+1) - 64*j)) + (1ULL<<((16+1)*(4*j+2) - 64*j)) + (1ULL<<((16+1)*(4*j+3) - 64*j)));
        _a[j] >>= (16-1);
        _a[j] &= (1ULL<<(16*0)) | (1ULL<<(16*1)) | (1ULL<<(16*2)) | (1ULL<<(16*3));
        _a[j] *= (1<<16) - 1;
    }

    for (i = 0; i < len; i++) {
        uint64_t b = (READ_INPUT(i) * ((1ULL<<(16*0)) + (1ULL<<(16*1)) + (1ULL<<(16*2)) + (1ULL<<(16*3))));

        #pragma GCC unroll 4
        for (j = 0; j < 4; j++) {
            perf_col[i][j] ^= _a[j] & b;
        }
    }
#elif defined(FEC_HAS_32BIT)
    int j;
    uint32_t _a[8];
    
    #pragma GCC unroll 8
    for (j = 0; j < 8; j++) {
        _a[j] = ((uint16_t)(((int16_t)(a << (16 - 1 - (j*2 + 0)))) >> 15)) | (((uint32_t)((uint16_t)(((int16_t)(a << (16 - 1 - (j*2 + 1)))) >> 15))) << 16);
    }

    for (i = 0; i < len; i++) {
        uint16_t b = READ_INPUT(i);
        uint32_t b_2 = b | (b << 16);

        #pragma GCC unroll 8
        for (j = 0; j < 8; j++) {
            perf_col[i][j] ^= _a[j] & b_2;
        }
    }
#endif
}

static inline void NORM_FUNC_NAME(const fec_perf_int_t* restrict perf_col, LEN_PARAM_TYPE len, OUTPUT_ARGS) {
    LEN_PARAM_TYPE i;

#if defined(FEC_HAS_CLMUL64)
    _poly_t _poly = _POLY_1VAL(POLY_G);

    for (i = 0; i < len - 1; i += 2) {

        _poly_t _c = _POLY_2VAL(perf_col[i], perf_col[i+1]);
        _poly_t _d = _POLY_VEC_SHIFT(_c, 16);
        _d = _POLY_CLMUL(_d, _poly);
        _c ^= _d;
        _d = _POLY_VEC_SHIFT(_d, 16);
        _d = _POLY_CLMUL(_d, _poly);
        _c ^= _d;
        
        WRITE_OUTPUT(i, _POLY_EXTRACT(_c, uint16_t, 0));
        WRITE_OUTPUT(i+1, _POLY_EXTRACT(_c, uint16_t, 2));
    }

    if (i == len - 1) {
        _poly_t _c = _POLY_1VAL(perf_col[i]);
        _poly_t _d = _c >> 16;
        _d = _POLY_CLMUL(_d, _poly);
        _c ^= _d;
        _d >>= 16;
        _d = _POLY_CLMUL(_d, _poly);
        _c ^= _d;
        
        WRITE_OUTPUT(i, _POLY_EXTRACT(_c, uint16_t, 0));
    }
#elif defined(FEC_HAS_CLMUL32)
    _poly_t _poly = _POLY_1VAL(POLY_G);

    for (i = 0; i < len; i++) {
        _poly_t _c = _POLY_1VAL(perf_col[i]);
        _poly_t _d = _c >> 16;
        _d = _POLY_CLMUL(_d, _poly);
        _c ^= _d;
        _d >>= 16;
        _d = _POLY_CLMUL(_d, _poly);
        _c ^= _d;
        
        WRITE_OUTPUT(i, _POLY_EXTRACT(_c, uint16_t, 0));
    }
#elif defined(FEC_HAS_128_INT_VEC)
    u32x8 shifts1 = {0, 2, 4, 6, 8, 10, 12, 14};
    u32x8 shifts2 = {1, 3, 5, 7, 9, 11, 13, 15};

    for (i = 0; i < len; i++) {
        u16x16 res = perf_col[i];
        u32x8 res1 = (((u32x8)res) & 0xFFFF) << shifts1;
        u32x8 res2 = (((u32x8)res) >> 16) << shifts2;

        u32x8 res3 = res1 ^ res2;
        u32x4 res4 = (u32x4)(my_mm256_extracti128_si256((__m256i)res3, 0) ^ my_mm256_extracti128_si256((__m256i)res3, 1));
        u32x4 res5 = res4 ^ my_mm_shuffle_epi32(res4, 1, 2, 3, 0);
        u32x4 res6 = res5 ^ my_mm_shuffle_epi32(res5, 2, 3, 0, 1); // res6 holds in each cell the int32 res

        u32x4 res7 = (res6 >> 16) << ((u32x4){0, 1, 3, 5}); // mul upper by POLY_G
        
        u32x4 res8 = res7 ^ my_mm_shuffle_epi32(res7, 1, 2, 3, 0);
        u32x4 res9 = res8 ^ my_mm_shuffle_epi32(res8, 2, 3, 0, 1);

        u32x4 res10 = (res9 >> 16) << ((u32x4){0, 1, 3, 5}); // mul upper by POLY_G
        
        u32x4 res11 = res10 ^ my_mm_shuffle_epi32(res10, 1, 2, 3, 0);
        u32x4 res12 = res11 ^ my_mm_shuffle_epi32(res11, 2, 3, 0, 1);

        WRITE_OUTPUT(i, ((u16x8)(res6 ^ res9 ^ res12))[0]);

    }
#else
    int j;
    for (i = 0; i < len; i++) {

        uint32_t res = 0;
        for (j = 0; j < 16; j++) {
            res ^= ((uint32_t)(((const uint16_t*)perf_col)[i*16 + j])) <<
#if defined(FEC_HAS_64_INT_VEC) || !defined(FEC_HAS_64BIT)
            j
#else
            (16 - 1 - j)
#endif
            ;
        }

        uint16_t res4 = res;
        uint32_t carry = res >> 16;

        carry ^= (carry << 1) ^ (carry << 3) ^ (carry << 5);

        res4 ^= carry;

        carry = carry >> 16;
        
        res4 ^= carry ^ (carry << 1) ^ (carry << 3) ^ (carry << 5);

        WRITE_OUTPUT(i, res4);
    }
#endif
}
