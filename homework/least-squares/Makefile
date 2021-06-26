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
	set ylabel "log(activity) a.u.";\
	set tics out;\
	set key out;\
	set title "least-squares fit F_c(x)=c_0+c_1*t to log(Rutherford/Soddy data)" ;\
	plot "$<" index 0 with errorbars title "data"\
	,"$<" index 1 using 1:2 with lines title "F_c(x)"\
	,"$<" index 1 using 1:3 with lines title "F_C(x)_{c0+dc0}"\
	,"$<" index 1 using 1:4 with lines title "F_C(x)_{c0-dc0}"\
	,"$<" index 2 using 1:3 with lines title "F_C(x)_{c1+dc1}"\
	,"$<" index 2 using 1:4 with lines title "F_C(x)_{c1-dc1}"\
	'|gnuplot


data.txt: main
	./$< > out.txt

main: main.o decomp.o

main.o decomp.o: decomp.h

clean:
	$(RM) *.o main out.txt data.txt *.svg
