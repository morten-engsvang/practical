CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: plot.svg

plot.svg: data.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "Time t (days)";\
	set ylabel "log(activity) arb. u.";\
	set tics out;\
	set key out;\
	set title "least-squares fit F_c(x)=c_0+c_1*t to log(Rutherford/Soddy data)" ;\
	plot "$<" index 0 with errorbars title "data"\
	,"$<" index 1 using 1:2 with lines title "F_c(x)"\
	'|gnuplot


data.txt: main
	./$< > data.txt

main: main.o decomp.o

main.o decomp.o: decomp.h

clean:
	$(RM) *.o main data.txt *.svg
