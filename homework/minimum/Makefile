CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt figure.svg

figure.svg: plotdata.txt Makefile
	echo '\
	set terminal svg background rgb "white";;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "Energy, E [GeV]";\
	set ylabel "Cross section, s [unspecified unit]";\
	set title "Higgs Boson data";\
	plot \
	 "$<" index 0 using 1:2:3 with yerrorbars title "Data"\
	,"$<" index 1 using 1:2 with line lt 2 title "Breit-Wigner fit"\
	'| gnuplot

plotdata.txt: main

out.txt: main
	./$< > out.txt

main: 


clean:
	$(RM) *.o main out.txt figure.svg plotdata.txt
