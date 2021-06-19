CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)
length = $(shell wc -l input.txt | cut -f1 -d ' ')
z = 4
.PHONY: plot.png

plot.png: default
	#awk '{print $$1, $$2}' out_data.txt > line1
	#awk '{print $$1, $$2}' input.txt > line2
	graph \
	--output-format png --bitmap-size 1920x1080 \
        --x-label "x" --y-label "y" \
        --top-label "Input points and interpolated points" \
	--line-mode 0 out_data.txt \
        --line-mode 0 --toggle-use-color --symbol 3 input.txt \
	> $@
	$(RM) line

default: main
	./$< input.txt $(length) 4

main: main.o binsearch.o

main.o binsearch.o: binsearch.h


clean:
	$(RM) *.o main out_data.txt out.txt *.png
