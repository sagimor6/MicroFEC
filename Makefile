
CC = gcc

CFLAGS += -pedantic-errors -Wno-pedantic -Wall -Wextra -Wnull-dereference -Waggressive-loop-optimizations -Wstrict-overflow=5 -Wuninitialized -fsanitize=undefined -fsanitize=address -g
CFLAGS += -DFEC_MIN_MEM

.PHONY all:
all: fec_test

.PHONY fec_test:
fec_test: fec_test.elf
	./fec_test.elf

fec_test.elf: main.c my_fec.c
	$(CC) $(CFLAGS) -o $@ $^

.PHONY clean:
clean:
	rm -f fec_test.elf
