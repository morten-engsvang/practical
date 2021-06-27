CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)


default: plot_plain.svg plot_halton.svg plot_plain2.svg plot_halton2.svg

plot_plain.svg: plot.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "N";\
	set ylabel "Error";\
	set tics out;\
	set key top right;\
	set yrange [0:0.4];\
	set xrange [0:500];\
	set title "Integral 1: Error vs. amount of sampling points N, x and y-axis has been truncated" ;\
	plot "$<" using 1:2 with lines title "plainmc estimated error"\
	,"$<" using 1:3 with lines title "plainmc actual error"\
	,"fit.txt" using 1:2 with lines title "Fit of 1/sqrt(N) to plainmc estimated error"\
	,"fit.txt" using 1:3 with lines title "Fit of 1/sqrt(N) to plainmc actual error"\
	'|gnuplot
	
plot_plain2.svg: plot2.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "N";\
	set ylabel "Error";\
	set tics out;\
	set key top right;\
	set title "Integral 2: Error vs. amount of sampling points N" ;\
	plot "$<" using 1:2 with lines title "plainmc estimated error"\
	,"$<" using 1:3 with lines title "plainmc actual error"\
	,"fit2.txt" using 1:2 with lines title "Fit of 1/sqrt(N) to plainmc estimated error"\
	,"fit2.txt" using 1:3 with lines title "Fit of 1/sqrt(N) to plainmc actual error"\
	'|gnuplot

plot_halton.svg: plot.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "N";\
	set ylabel "Error";\
	set tics out;\
	set key top right;\
	set yrange [0:0.4];\
	set title "Integral 1: Error vs. amount of sampling points N, y-axis has been truncated" ;\
	plot "$<" using 1:4 with lines title "haltonmc estimated error"\
	,"$<" using 1:5 with lines title "haltonmc actual error"\
	,"fit.txt" using 1:4 with lines title "Fit of 1/sqrt(N) to haltonmc estimated error"\
	,"fit.txt" using 1:5 with lines title "Fit of 1/sqrt(N) to haltonmc actual error"\
	'|gnuplot

plot_halton2.svg: plot2.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "N";\
	set ylabel "Error";\
	set tics out;\
	set key top right;\
	set title "Integral 2: Error vs. amount of sampling points N" ;\
	plot "$<" using 1:4 with lines title "haltonmc estimated error"\
	,"$<" using 1:5 with lines title "haltonmc actual error"\
	,"fit2.txt" using 1:4 with lines title "Fit of 1/sqrt(N) to haltonmc estimated error"\
	,"fit2.txt" using 1:5 with lines title "Fit of 1/sqrt(N) to haltonmc actual error"\
	'|gnuplot

plot.txt: main
	./$< > out.txt

main: main.o montecarlo.o minimum.o

main.o montecarlo.o minimum.o: montecarlo.h minimum.h

clean:
	$(RM) *.o main out.txt fit.txt fit2.txt plot.txt plot2.txt *.err *.svg
