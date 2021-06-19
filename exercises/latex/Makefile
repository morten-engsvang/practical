CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

exp.pdf: exp.tex graph.png
	pdflatex $<
	pdflatex $<

graph.png: data.txt Makefile
	awk '{print $$1, $$2}' $< > line1
	awk '{print $$1, $$3}' $< > line2
	graph \
	--output-format png --bitmap-size 1024x768 \
        --x-label "x" --y-label "y" \
        --top-label "Exponential functions" \
        --line-mode 1 line1 \
        --line-mode 2 line2 \
	> $@
	$(RM) line*

data.txt: main
	./$< > $@


clean: 
	$(RM) exp.pdf *.log *.aux *.o main *.png
