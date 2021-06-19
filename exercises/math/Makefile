CC = gcc

default: out.txt
	cat out.txt

out.txt: math
	./math > out.txt

math: math.o
	$(CC) -o math math.o

math.o: math.c
	$(CC) -c math.c

clean:
	rm -f math.o math out.txt
