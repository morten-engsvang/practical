CC = gcc
CFLAGS = -std=gnull
LDLIBS = -lm


default: out.txt
	cat out.txt

out.txt: epsilon
	./epsilon > out.txt

epsilon: epsilon.o equal.o
	$(CC) -o epsilon epsilon.o equal.o $(LDLIBS)

epsilon.o: epsilon.c
	$(CC) -c epsilon.c -fwrapv -O

equal.o: equal.c
	$(CC) -c equal.c

clean:
	rm -f epsilon.o separat.o out.txt epsilon
