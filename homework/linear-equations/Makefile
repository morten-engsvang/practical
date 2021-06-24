CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

time.svg: default
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "n";\
	set ylabel "time [seconds]";\
	set tics out;\
	set key top right;\
	set title "Time vs. matrix size" ;\
	plot "outtime.txt" using 1:2 with lines title "My time"\
	,"outtime.txt" using 1:3 with lines title "gsl linalg QR decomp"\
	'|gnuplot

default: main
	./$< > out.txt

main:

clean:
	$(RM) *.o main out.txt

