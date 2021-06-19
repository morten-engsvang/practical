CLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)
.PHONY: out.txt

out.txt: main
	./$< > $@

main: 

clean:
	$(RM) *.o main *.txt
