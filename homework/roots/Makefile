CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.svg

out.svg: numode.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "r";\
	set ylabel "y";\
	set tics out;\
	set key right top;\
	set title "Radial function for the S-state from r = 0.0001 to r = 8" ;\
	plot "$<" using 1:2 with lines title "Numerical solution, epsilon = -0.500025"\
	,"exactode.txt" using 1:2 with lines title "Exact solution, epsilon = -1/2"\
	'|gnuplot


numode.txt: main
	./$< > out.txt


main: main.o qr.o rkstep12.o

main.o qr.o rkstep12.o: qr.h rkstep12.h

clean:
	$(RM) *.o main *.txt *.png *.svg
