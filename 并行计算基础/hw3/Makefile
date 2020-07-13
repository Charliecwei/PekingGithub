CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = hw3
all: ${EXE}


Hw3: omp_hw3ss.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm


clean:
	rm -f $(OBJ) $(EXE)

