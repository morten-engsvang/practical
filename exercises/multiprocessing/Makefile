CFLAGS = -Wall -O -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp
.PHONY: default

default: main main_omp

main:

main_omp:

clean:
	$(RM) *.o main main_omp
