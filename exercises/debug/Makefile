CFLAGS = -o -w -std=gnu99
CLAGS += `shell gsl-config --cflags`
LDLIBS += $(shell gsl-config --libs)

all: output ; cat $<
output: main ; ./$< > $@
clean: ; $(RM) main output
