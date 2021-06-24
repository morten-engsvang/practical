CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)


default: plot.svg
	
plot.svg: neural.txt function.txt params.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set tics out;\
	set key right bottom;\
	set title "Neural network (5 neurons) approximating cos(5x-1)*exp(-x*x)" ;\
	plot "neural.txt" using 1:2 with lines title "Neural network approximation"\
	,"function.txt" using 1:2 with points ls 6 title "Analytical solution"\
	,"neural_dif.txt" using 1:2 with lines title "Neural network derivative"\
	,"neural_int.txt" using 1:2 with lines title "Neural network antiderivative"\
	'|gnuplot



neural.txt function.txt params.txt: main
	./$< > out.txt

main: main.o qnewton.o ann.o

main.o qnewton.o ann.o: ann.h qnewton.h

clean:
	$(RM) *.o main *.txt *.svg
