CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: main
	./$< > out.txt

main:

clean:
	$(RM) *.o main out.txt

