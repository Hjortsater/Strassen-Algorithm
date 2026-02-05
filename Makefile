CC := gcc

HAVE_CLANG_OMP := $(shell \
	echo 'int main(){}' | \
	clang -fopenmp -lomp -x c - -o /dev/null 2>/dev/null && \
	echo yes)

ifeq ($(HAVE_CLANG_OMP),yes)
    CC := clang
    OMP_CFLAGS := -fopenmp
    OMP_LDFLAGS := -lomp
    $(info Using clang + libomp)
else
    CC := gcc
    OMP_CFLAGS := -fopenmp
    OMP_LDFLAGS :=
    $(info Using gcc + OpenMP)
endif

CFLAGS = -Wall -Wextra -O3 \
         $(OMP_CFLAGS) \
         -DPRINT=0 \
         -DOPTMULT=1 \
         -DOPTCUTOFF=32 \
         -DPARALLEL_CUTOFF=64

SRC = main.c matrix_alg.c strassen.c
TARGET = main

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(OMP_LDFLAGS) -o $(TARGET)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)
