CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: main
	./$< > out.txt

main: main.o montecarlo.o

main.o montecarlo.o: montecarlo.h

clean:
	$(RM) *.o main *.txt *.err 
