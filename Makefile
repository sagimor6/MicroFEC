
#CROSS_PREFIX = /home/sagi/git_repos/static-python/aarch64--musl--stable-2022.08-2/bin/aarch64-buildroot-linux-musl-
#CROSS_PREFIX = /home/sagi/git_repos/static-python/aa/static-python/armv7-eabihf--musl--bleeding-edge-2023.11-1/bin/arm-buildroot-linux-musleabihf-
#CROSS_RUNNER = qemu-aarch64-static
#CROSS_RUNNER = qemu-arm-static
CC = $(CROSS_PREFIX)clang
OBJCOPY = $(CROSS_PREFIX)objcopy

CFLAGS += -pedantic-errors -Wno-pedantic -Wall -Wextra -Wnull-dereference -Waggressive-loop-optimizations -Wstrict-overflow=5 -Wuninitialized
#CFLAGS += -DFEC_MIN_MEM -DFEC_LARGE_K

DEBUG_CFLAGS += -g
#DEBUG_CFLAGS += -fsanitize=undefined -fsanitize=address 

OPT_CFLAGS += -flto -ffat-lto-objects -ffunction-sections -fdata-sections -O2 -fvisibility=hidden
OPT_CFLAGS += -mpclmul
OPT_CFLAGS += -mavx2 -mtune=cannonlake
# OPT_CFLAGS += -mavx
# OPT_CFLAGS += -m32 -mno-sse2 -mno-sse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-avx2 -mno-avx -mno-sse -mno-mmx
# OPT_CFLAGS += -march=armv8-a+crypto
# OPT_CFLAGS += -march=armv8-a+crypto -mfpu=crypto-neon-fp-armv8
OPT_LDFLAGS = -Wl,--gc-sections

#CROSS_FLAGS += -static
#VALGRIND = valgrind -s --show-leak-kinds=all --partial-loads-ok=no --expensive-definedness-checks=yes

BUILD = build

.PHONY all:
all: $(BUILD)/libmicro_fec.so $(BUILD)/libmicro_fec.a fec_test

.PHONY fec_test:
fec_test: $(BUILD)/fec_test.elf
	$(VALGRIND) $(CROSS_RUNNER) $(BUILD)/fec_test.elf 1000 100 1400

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
