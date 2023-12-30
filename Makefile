
CC = gcc
OBJCOPY = objcopy

CFLAGS += -pedantic-errors -Wno-pedantic -Wall -Wextra -Wnull-dereference -Waggressive-loop-optimizations -Wstrict-overflow=5 -Wuninitialized
CFLAGS += -DFEC_MIN_MEM -DFEC_LARGE_K

DEBUG_CFLAGS += -g
#DEBUG_CFLAGS += -fsanitize=undefined -fsanitize=address 

OPT_CFLAGS += -flto -ffat-lto-objects -ffunction-sections -fdata-sections -Os -fvisibility=hidden -s
OPT_CFLAGS += -mpclmul
OPT_LDFLAGS = -Wl,--gc-sections

BUILD = build

.PHONY all:
all: $(BUILD)/libmicro_fec.so $(BUILD)/libmicro_fec.a fec_test

.PHONY fec_test:
fec_test: $(BUILD)/fec_test.elf
	$(BUILD)/fec_test.elf 1000 100 1400

$(BUILD)/:
	mkdir -p $@

$(BUILD)/fec_test.elf: main.c micro_fec.c | $(BUILD)/
	$(CC) $(CFLAGS) $(DEBUG_CFLAGS) $(OPT_CFLAGS) $(OPT_LDFLAGS) -o $@ $^

$(BUILD)/micro_fec.o_static: micro_fec.c | $(BUILD)/
	$(CC) $(CFLAGS) $(OPT_CFLAGS) -c -o $@ $^

$(BUILD)/libmicro_fec.so: micro_fec.c | $(BUILD)/
	$(CC) $(CFLAGS) $(OPT_CFLAGS) $(OPT_LDFLAGS) -D_FEC_DO_EXPORTS -shared -fpic -o $@ $^

$(BUILD)/libmicro_fec.a: $(BUILD)/micro_fec.o_static | $(BUILD)/
	rm -f $@
	ar rcs $@ $^

.PHONY clean:
clean:
	[ ! -d $(BUILD) ] || rm -r $(BUILD)
