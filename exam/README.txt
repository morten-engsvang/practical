This set of files has been made by Morten Engsvang: Student Number: 201808926
26%22 = 4 which means that I had been assigned exam project 4.

PROBLEM STATEMENT:
- 4 - Multidimensional pseudo-random (plain Monte Carlo) vs. quasi-random (Halton and/or lattice sequence) integrators.

Investigate the convergence rates (of some interesting integrals in different dimensions) as a function of the number of sample points.

----------

When cleaned (make clean) this folder contains the following:
	main.c	:	Main driving programs produces: out.txt, plot.txt, plot2.txt, fit.txt and fit2.txt
	minimum.c:	Contains the quasi-newton minimisation algorith from the minimisation homework
	minimum.h:	The corresponding header file to minimum.c
	montecarlo.c:	Contains implementations of the plain Monte Carlo algorithm and the quasi-random Monte Carlo Algorithm using the Halton Sequence.
	Makefile:	Creates plot_plain.svg, plot_plain2.svg, plot_halton.svg and plot_halton2.svg
	README.txt:	The file currently being read :)

Comments on the project:
	The approach which I took was to calculate the integrals using the already implemented Monte-Carlo algorithms. 
	In these I had implemented estimates of the error which I used alongside the actual error that could be calculated because the real value is known.
	I then used the minimisation algorithm previously implemented in a homework in order to fit f(x)=a+b*1/sqrt(x)
	because it is known that the error from a random walk scales as 1/sqrt(x). This is the purpose of the defined errorfunctions.
	This was only done for two integrals:
	Integral 1: The 'difficult and singular' integral from the Monte Carlo Integration homework part A.
	Integral 2: 4*sqrt(1-x^2) from 0 to 1 from the Adaptive Integration homework part A.
	I did not do it for more integral because I do not know any more interesting integrals.

Comments on the code:
	The structure of the code follows the overall structure of the project. Therefore I will not comment further on that.
	What I feel is missing in the code:
	1: The code is not very modular i.e. it is hard to implement evaluation of other integrals.
	2: I evaluate very few integrals, however the amount of dimensions are different between the two that was evaluated.
	3: I could also have implemented the Monte Carlo algorithm using either the lattice sequence or one using recursive stratified sampling.
	4: Prettier code in general.

Comments on the results:
	The results can be viewed in out.txt and the .svg files with more details to be found.
	It can be seen that the 1. integral, which is singular, does not converge very nicely which also makes the fit relatively bad.
	The 2. integral behaves much more nicely and it can be seen that the fit describes the trend in convergence relatively well.
	It can also be seen that the change in sampling between plain and halton Monte Carlo makes the convergence much more
	
Self-evaluation:
	I feel like I have technically fulfilled the project describtion, however due to; the low amount of integrals studied,
	only looking at two sampling techniques, the code quality, and only attempting to fit one function, I feel that the project is lacking some work.
	Therefore I would give this around 5-6 out of 10. 
	Or in homework terms: part A has been done.
	
	
