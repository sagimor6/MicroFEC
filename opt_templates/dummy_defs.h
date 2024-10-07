

#ifndef FMA_FUNC_NAME
#define FMA_FUNC_NAME fec_col_fma
#endif

#ifndef NORM_FUNC_NAME
#define NORM_FUNC_NAME fec_norm
#endif

#ifndef INIT_FUNC_NAME
#define INIT_FUNC_NAME fec_init_col
#ifndef INIT_FUNC_ARGS
#define INIT_FUNC_ARGS fec_int_t *__dummy
#endif

#ifndef INIT_FUNC_READ_INPUT
#define INIT_FUNC_READ_INPUT(j) __dummy[j]
#endif
#endif

#ifndef LEN_PARAM_TYPE
#define LEN_PARAM_TYPE size_t
#endif

#ifndef INPUT_ARGS
#define INPUT_ARGS const fec_int_t *__dummy
#endif

#ifndef OUTPUT_ARGS
#define OUTPUT_ARGS fec_int_t *__dummy
#endif

#ifndef READ_INPUT
#define READ_INPUT(j) __dummy[j]
#endif

#ifndef WRITE_OUTPUT
#define WRITE_OUTPUT(j, val) __dummy[j] = val
#endif


