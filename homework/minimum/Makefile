CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: main

main: main.o qr.o rkstep12.o

main.o qr.o rkstep12.o: qr.h rkstep12.h

clean:
	$(RM) *.o main *.txt *.png *.svg
