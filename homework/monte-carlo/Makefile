CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: plot.svg plot_actual.svg

plot.svg: plot.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "N";\
	set ylabel "Error";\
	set tics out;\
	set key top right;\
	set title "Error vs. amount of sampling points N" ;\
	plot "$<" using 1:2 with lines title "plainmc estimated error"\
	,"$<" using 1:4 with lines title "haltonmc estimated error"\
	'|gnuplot

plot_actual.svg: plot.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "N";\
	set ylabel "Error";\
	set tics out;\
	set key top right;\
	set title "Error vs. amount of sampling points N" ;\
	plot "$<" using 1:3 with lines title "plainmc actual error"\
	,"$<" using 1:5 with lines title "haltonmc actual error"\
	'|gnuplot

plot.txt: main
	./$< > out.txt

main: main.o montecarlo.o

main.o montecarlo.o: montecarlo.h

clean:
	$(RM) *.o main *.txt *.err *.svg
