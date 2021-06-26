CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)


default: eigen.png energy.png

#Plotting in the same way as my mate Kasper Larsen
energy.png: energy.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$n$$";\
	set ylabel "$$E(n)$$ in units of $$2mL^2/\hbar^2$$";\
	set xrange [0:16];\
	set title "Quantum Particle in a Box";\
	plot \
	 "$<" using 1:2 with points pt 2 ps 0.3 title "calculated values"\
	,"$<" using 1:3 with points pt 1 ps 0.3 title "exact values"\
	'| pyxplot


eigen.png: eigen.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$\xi$$";\
	set ylabel "$$u_n(\xi)$$";\
	set title "Eigenfunctions";\
	plot \
	 "$<" index 0 using 1:2 with line lt 2 color blue title "$$n=1$$ calculated eigenfunction"\
	,"$<" index 0 using 1:3 with line lt 1 lw 0.5 color blue title "$$n=1$$ exact eigenfunction"\
	,"$<" index 1 using 1:2 with line lt 2 color red title "$$n=2$$ calculated eigenfunction"\
	,"$<" index 1 using 1:3 with line lt 1 lw 0.5 color red title "$$n=2$$ exact eigenfunction"\
	,"$<" index 2 using 1:2 with line lt 2 color black title "$$n=3$$ calculated eigenfunction"\
	,"$<" index 2 using 1:3 with line lt 1 lw 0.5 color black title "$$n=3$$ exact eigenfunction"\
	'| pyxplot

eigen.txt: main 
	./$< > out.txt

main:

clean:
	$(RM) *.o main *.txt *.png *.svg
