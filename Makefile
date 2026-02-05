CC = clang

# aggressive optimization for speed: vectorize, tune, link-time-opt, fast-math
CFLAGS = -Wall -Wextra -Ofast -march=native -flto -funroll-loops -fopenmp -ffast-math \
         -DPRINT=0 \
         -DOPTMULT=2 \
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
