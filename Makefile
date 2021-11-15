CFLAGS := ${CFLAGS} -std=c99 -g

all:
	$(CC) $(CFLAGS) -c quadtree.c
	$(CC) $(CFLAGS) -c spec.c
	$(CC) $(CFLAGS) butterfly.c -o butterfly quadtree.o spec.o

	$(CC) $(CFLAGS) build_quadtree.c -o build_quadtree quadtree.o
