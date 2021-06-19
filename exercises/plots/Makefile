CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS = $(shell gsl-config --libs)
.PHONY: default

default: erf_graph.png gamma_graph.png

erf_graph.png: erf_data.txt Makefile 
	awk '{print $$1, $$2}' $< > line1
	awk '{print $$1, $$3}' $< > line2
	awk '{print $$1, $$4}' $< > line3
	graph \
	--output-format png --bitmap-size 1024x768 \
	--x-label "x" --y-label "y" \
	--top-label "Error function" \
	--line-mode 1 line1 \
	--line-mode 2 line2 \
	--line-mode 3 line3 \
	> $@
	$(RM) line*

gamma_graph.png: gamma_data.txt Makefile
	awk '{print $$1, $$2}' $< > line1
	awk '{print $$1, $$3}' $< > line2
	awk '{print $$1, $$4}' $< > line3
	graph \
	--output-format png --bitmap-size 1024X768 \
	--x-label "x" --y-label "y" \
	--top-label "Gamma function" \
	--line-mode 1 line1 \
	--line-mode 2 line2 \
	--line-mode 3 line3 \
	> $@
	$(RM) line*

gamma_data.txt: gamma
	./$< > $@

erf_data.txt: erf
	./$< > $@

gamma: gamma.o mygamma.o

gamma.o mygamma.o: mygamma.h

erf: erf.o myerf.o

erf.o myerf.o: myerf.h

clean:
	$(RM) *.o erf gamma erf_data.txt gamma_data.txt *.png
