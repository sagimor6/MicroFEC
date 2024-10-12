
#CROSS_PREFIX = /home/sagi/git_repos/static-python/aarch64--musl--stable-2022.08-2/bin/aarch64-buildroot-linux-musl-
#CROSS_PREFIX = /home/sagi/git_repos/static-python/aa/static-python/armv7-eabihf--musl--bleeding-edge-2023.11-1/bin/arm-buildroot-linux-musleabihf-
#CROSS_RUNNER = qemu-aarch64-static
#CROSS_RUNNER = qemu-arm-static
CC_NAME ?= clang
CC = $(CROSS_PREFIX)$(CC_NAME)
OBJCOPY = $(CROSS_PREFIX)objcopy

check_cc_flag = $(if $(shell if $(CC) $(1) -Werror -Wl,--help -xc -o /dev/null /dev/null >/dev/null 2>/dev/null; then echo a; fi; ),$(1),)

CFLAGS += -pedantic-errors -Wno-pedantic -Wall -Wextra -Wnull-dereference -Wstrict-overflow=5 -Wuninitialized
CFLAGS += -Waggressive-loop-optimization
#CFLAGS += -DFEC_MIN_MEM -DFEC_LARGE_K

DEBUG_CFLAGS += -g
#DEBUG_CFLAGS += -fsanitize=undefined -fsanitize=address 

OPT_CFLAGS += -ffunction-sections -fdata-sections -O3 -fvisibility=hidden
OPT_CFLAGS += -flto -fuse-linker-plugin -ffat-lto-objects -flto=auto -flto-partition=one
OPT_CFLAGS += -falign-loops=32
OPT_LDFLAGS += -Wl,--gc-sections

#ARCH_CFLAGS += -march=native

check_option = $(shell echo | $(CC) $(ARCH_CFLAGS) -v -E - 2>&1 >/dev/null | grep -E -e '\s$(1)\s')
has_jcc_erratum := $(or $(has_jcc_erratum),$(call check_option,-march=skylake),$(call check_option,-march=skylake-avx512),$(call check_option,-march=cascadelake))
has_jcc_erratum := $(or $(has_jcc_erratum),$(call check_option,-mtune=skylake),$(call check_option,-mtune=skylake-avx512),$(call check_option,-mtune=cascadelake))
has_jcc_erratum := $(or $(has_jcc_erratum),$(call check_option,-target-cpu skylake),$(call check_option,-target-cpu skylake-avx512),$(call check_option,-target-cpu cascadelake))
has_jcc_erratum := $(or $(has_jcc_erratum),$(call check_option,-tune-cpu skylake),$(call check_option,-tune-cpu skylake-avx512),$(call check_option,-tune-cpu cascadelake))
jcc_erratum_cflags = -mbranches-within-32B-boundaries -Wa,-mbranches-within-32B-boundaries
ifneq ($(has_jcc_erratum),)
ARCH_CFLAGS += $(foreach cflag,$(jcc_erratum_cflags),$(call check_cc_flag,$(cflag)))
endif

CFLAGS += -DPERF_DEBUG

# 32bit:
# ARCH_CFLAGS += -m32

# pclmul:
# ARCH_CFLAGS += -mpclmul
# avx2:
# ARCH_CFLAGS += -mno-pclmul -mavx2
# avx:
# ARCH_CFLAGS += -mno-pclmul -mno-avx2 -mavx
# sse2:
# ARCH_CFLAGS += -mno-pclmul -mno-avx2 -mno-avx -msse2
# sse:
# ARCH_CFLAGS += -mno-pclmul -mno-avx2 -mno-avx -mno-sse4.2 -mno-sse4.1 -mno-sse4 -mno-sse3 -mno-sse2 -msse
# mmx:
# ARCH_CFLAGS += -mno-pclmul -mno-avx2 -mno-avx -mno-sse4.2 -mno-sse4.1 -mno-sse4 -mno-sse3 -mno-sse2 -mno-sse -mmmx
# 64bit:
# ARCH_CFLAGS += -mno-pclmul -mno-avx2 -mno-avx -mno-sse4.2 -mno-sse4.1 -mno-sse4 -mno-sse3 -mno-sse2 -mno-sse -mno-mmx
# 32bit:
# ARCH_CFLAGS += -mno-pclmul -mno-avx2 -mno-avx -mno-sse4.2 -mno-sse4.1 -mno-sse4 -mno-sse3 -mno-sse2 -mno-sse -mno-mmx -m32

# aarch64 clmul:
# ARCH_CFLAGS += -march=armv8-a+crypto

# arm clmul:
# ARCH_CFLAGS += -march=armv8-a+crypto -mfpu=crypto-neon-fp-armv8



#CROSS_FLAGS += -static
#VALGRIND = valgrind -s --show-leak-kinds=all --partial-loads-ok=no --expensive-definedness-checks=yes

OPT_CFLAGS := $(foreach cflag,$(OPT_CFLAGS),$(call check_cc_flag,$(cflag)))
CFLAGS := $(foreach cflag,$(CFLAGS),$(call check_cc_flag,$(cflag)))
OPT_LDFLAGS := $(foreach cflag,$(OPT_LDFLAGS),$(call check_cc_flag,$(cflag)))

CFLAGS += $(ARCH_CFLAGS)

TEST_PARAMS ?= 10000 2000 500

BUILD = build

.PHONY all:
all: $(BUILD)/libmicro_fec.so $(BUILD)/libmicro_fec.a fec_test

.PHONY fec_test:
fec_test: $(BUILD)/fec_test.elf
	sudo $(VALGRIND) $(CROSS_RUNNER) $(BUILD)/fec_test.elf $(TEST_PARAMS)

$(BUILD)/:
	mkdir -p $@

$(BUILD)/fec_test.elf: main.c micro_fec.c | $(BUILD)/
	$(CC) $(CFLAGS) $(DEBUG_CFLAGS) $(OPT_CFLAGS) $(OPT_LDFLAGS) $(CROSS_FLAGS) -o $@ $^

$(BUILD)/micro_fec.o_static: micro_fec.c | $(BUILD)/
	#-flinker-output=nolto-rel
	$(CC) $(CFLAGS) $(OPT_CFLAGS) -Wl,--gc-sections -Wl,--gc-keep-exported -D_FEC_DO_EXPORTS -r -o $@ $^
	#$(CC) $(CFLAGS) $(OPT_CFLAGS) -Wl,--retain-symbols-file -Wl,exported_symbols.txt -r -o $@ $^
	#$(CC) $(CFLAGS) $(OPT_CFLAGS) -Wl,--version-script=verr -r -o $@ $^

# $(BUILD)/micro_fec.o_static2: $(BUILD)/micro_fec.o_static | $(BUILD)/
# 	$(OBJCOPY) --localize-hidden $^ $@

# $(BUILD)/micro_fec.o_static3: $(BUILD)/micro_fec.o_static2 | $(BUILD)/
# 	$(OBJCOPY) -x $^ $@

# $(BUILD)/micro_fec.o_static4: $(BUILD)/micro_fec.o_static3 | $(BUILD)/
# 	$(CC) $(CFLAGS) $(OPT_CFLAGS) -r -o $@ $^

$(BUILD)/libmicro_fec.so: micro_fec.c | $(BUILD)/
	$(CC) $(CFLAGS) $(OPT_CFLAGS) $(OPT_LDFLAGS) -D_FEC_DO_EXPORTS -shared -fpic -o $@ $^

$(BUILD)/libmicro_fec.a: $(BUILD)/micro_fec.o_static | $(BUILD)/
	rm -f $@
	ar rcs $@ $^

.PHONY clean:
clean:
	[ ! -d $(BUILD) ] || rm -r $(BUILD)
