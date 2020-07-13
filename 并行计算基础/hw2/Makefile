CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = cpi

all: ${EXE}

cpi: omp_hw2.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

clean:
	rm -f $(OBJ) $(EXE)

