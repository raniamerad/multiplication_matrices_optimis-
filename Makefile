CC=gcc
CFLAGS=-Wall -g

EXEC=matmul

.PHONY: clean

all: matmul

test: matmul
	./$(EXEC)

matmul: matmul.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(EXEC)
