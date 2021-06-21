CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: main
	./$< > out.txt


main: main.o jacobi.o

main.o jacobi.o: jacobi.h

clean:
	$(RM) *.o main *.txt *.png *.svg
