CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)
length = $(shell wc -l input.txt | cut -f1 -d ' ')
z = 4
default: plot.svg quad.svg int_afledt.svg

int_afledt.svg: input.txt Makefile 
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "t";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "Integral and derivative using interpolation" ;\
	plot "quadinteg.txt" using 1:2 with lines title "Quadratic integral"\
	,"quadderiv.txt" using 1:2 with lines title "Quadratic derivative"\
	,"lininteg.txt" using 1:2 with lines title "Linear Integral"\
	,"outlinteg.txt" using 1:2 with lines title "Linear Integral using GSL"\
	'|gnuplot


plot.svg: input.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "t";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "Input points and interpolated points" ;\
	plot "out_data.txt" using 1:2 with points title "Interpolation"\
	,"input.txt" using 1:2 with points ls 6 title "Input points"\
	,"outlspline.txt" using 1:2 with points title "GSL Interpolation"\
	'|gnuplot


quad.svg: input.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "t";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "Input points and quadratically interpolated points" ;\
	plot "outquadspline.txt" using 1:2 with lines title "Quadratic interpolation"\
	,"input.txt" using 1:2 with points ls 6 title "Input points"\
	'|gnuplot


input.txt: main
	./$< input.txt $(length) 4

main: main.o binsearch.o spline.o

main.o binsearch.o spline.o: binsearch.h


clean:
	$(RM) *.o main out_data.txt outlinteg.txt outlspline.txt out.txt *.png *.svg outquadspline.txt quadinteg.txt quadderiv.txt lininteg.txt
