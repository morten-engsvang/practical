CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS = $(shell gsl-config --libs)
.PHONY: default

default: integral_a erf_graph.png

erf_graph.png: b_data.txt Makefile
	awk '{print $$1, $$2}' $< > erf
	graph \
	--output-format png --bitmap-size 1024x768 \
	--x-label "x" --y-label "y" \
	--top-label "Error function via. integral" \
	--line-mode 1 erf \
	> $@
	$(RM) erf

b_data.txt: integral_b
	./$< > b_data.txt

integral_a:

integral_b:

clean:
	$(RM) *.o integral_a integral_b *.txt *.png
