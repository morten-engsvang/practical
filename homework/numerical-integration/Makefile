CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: main
	./$< > out.txt

main: main.o integrate.o

main.o integrate.o: integrate.h

clean:
	$(RM) *.o main *.txt *.err 
