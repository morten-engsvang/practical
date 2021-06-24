CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)
length = $(shell wc -l input.txt | cut -f1 -d ' ')
z = 4
.PHONY: plot.svg

plot.svg: default
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



default: main
	./$< input.txt $(length) 4

main: main.o binsearch.o

main.o binsearch.o: binsearch.h


clean:
	$(RM) *.o main out_data.txt outlinteg.txt outlspline.txt out.txt *.png *.svg
