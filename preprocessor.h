
#ifndef __PREPROCESSOR_H__
#define __PREPROCESSOR_H__

#define _GLUE(x,y) x##y
#define GLUE(x,y) _GLUE(x,y)

#define INC(x) _GLUE(INC_, x)
#define INC_0 1
#define INC_1 2
#define INC_2 3
#define INC_3 4
#define INC_4 5
#define INC_5 6
#define INC_6 7
#define INC_7 8
#define INC_8 9
#define INC_9 10
#define INC_10 11
#define INC_11 12
#define INC_12 13
#define INC_13 14
#define INC_14 15
#define INC_15 16
#define INC_16 17
#define INC_17 18
#define INC_18 19
#define INC_19 20
#define INC_20 21
#define INC_21 22
#define INC_22 23
#define INC_23 24
#define INC_24 25
#define INC_25 26
#define INC_26 27
#define INC_27 28
#define INC_28 29
#define INC_29 30
#define INC_30 31
#define INC_31 32

#define DEC(x) _GLUE(DEC_, x)
#define DEC_1 0
#define DEC_2 1
#define DEC_3 2
#define DEC_4 3
#define DEC_5 4
#define DEC_6 5
#define DEC_7 6
#define DEC_8 7
#define DEC_9 8
#define DEC_10 9
#define DEC_11 10
#define DEC_12 11
#define DEC_13 12
#define DEC_14 13
#define DEC_15 14
#define DEC_16 15
#define DEC_17 16
#define DEC_18 17
#define DEC_19 18
#define DEC_20 19
#define DEC_21 20
#define DEC_22 21
#define DEC_23 22
#define DEC_24 23
#define DEC_25 24
#define DEC_26 25
#define DEC_27 26
#define DEC_28 27
#define DEC_29 28
#define DEC_30 29
#define DEC_31 30
#define DEC_32 31

#define _ARG2(x, n, ...) n
#define _IS_PROBE(...) _ARG2(__VA_ARGS__, 0)
#define _PROBE() ~, 1

#define NOT(x) _IS_PROBE(_GLUE(_NOT_, x))
#define _NOT_0 _PROBE()

#define BOOL(x) NOT(NOT(x))

#define _IF(c) _GLUE(_IF_, c)
#define _IF_0(t, ...) __VA_ARGS__
#define _IF_1(t, ...) t

#define IF(c) _IF(BOOL(c))

#define EAT(...)
#define EXPAND(...) __VA_ARGS__
#define WHEN(c) IF(c)(EXPAND, EAT)

#define EMPTY()
#define DEFER(id) id EMPTY()
#define OBSTRUCT(id) id DEFER(EMPTY)()



#define EVAL(...) EVAL_512(__VA_ARGS__)
#define EVAL_512(...) EVAL_256(EVAL_256(__VA_ARGS__))
#define EVAL_256(...) EVAL_128(EVAL_128(__VA_ARGS__))
#define EVAL_128(...) EVAL_64(EVAL_64(__VA_ARGS__))
#define EVAL_64(...) EVAL_32(EVAL_32(__VA_ARGS__))
#define EVAL_32(...) EVAL_16(EVAL_16(__VA_ARGS__))
#define EVAL_16(...) EVAL_8(EVAL_8(__VA_ARGS__))
#define EVAL_8(...) EVAL_4(EVAL_4(__VA_ARGS__))
#define EVAL_4(...) EVAL_2(EVAL_2(__VA_ARGS__))
#define EVAL_2(...) EVAL_1(EVAL_1(__VA_ARGS__))
#define EVAL_1(...) __VA_ARGS__


#define REPEAT(count, macro, ...) \
    WHEN(count) \
    ( \
        OBSTRUCT(REPEAT_INDIRECT) () \
        ( \
            DEC(count), macro, __VA_ARGS__ \
        ) \
        OBSTRUCT(macro) \
        ( \
            DEC(count), __VA_ARGS__ \
        ) \
    )
#define REPEAT_INDIRECT() REPEAT

#endif // __PREPROCESSOR_H__
