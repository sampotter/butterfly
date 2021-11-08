CFLAGS := ${CFLAGS} -std=c99

all:
	$(CC) -c quadtree.c
	$(CC) -c spec.c
	$(CC) butterfly.c -o butterfly quadtree.o spec.o
