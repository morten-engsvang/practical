CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: test.svg sir.svg


test.svg: test.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "t";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "Plot of test function, u, with u(0)=0, du/dt=1, du^2/dt^2 = -u" ;\
	plot "$<" using 1:2 with lines title "y"\
	,"$<" using 1:3 with lines title "dy/dt"\
	'|gnuplot

sir.svg: sir.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "t";\
	set ylabel "amount";\
	set tics out;\
	set key out;\
	set title "The SIR model with TC = 3 and TR = 20" ;\
	plot "$<" using 1:2 with lines title "S"\
	,"$<" using 1:3 with lines title "I"\
	,"$<" using 1:4 with lines title "R"\
	'|gnuplot


test.txt: main
	./$<


main: main.o rkstep12.o

main.o rkstep12.o: rkstep12.h

clean:
	$(RM) *.o main *.txt *.png *.svg
