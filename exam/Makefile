CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)


default: main
	./$< > out.txt

main: main.o montecarlo.o least_squares.o

main.o montecarlo.o least_squares.o: montecarlo.h least_squares.h

clean:
	$(RM) *.o main *.txt *.err *.svg
