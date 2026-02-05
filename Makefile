CC = clang

CFLAGS = -Wall -Wextra -O3 -fopenmp\
         -DPRINT=0 \
         -DOPTMULT=1 \
         -DOPTCUTOFF=32 \
         -DPARALLEL_CUTOFF=64

SRC = main.c matrix_alg.c strassen.c

TARGET = main

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)
