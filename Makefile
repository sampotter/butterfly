CFLAGS := ${CFLAGS} -std=c99 -g -DBF_DEBUG -lm

all:
	$(CC) $(CFLAGS) -c quadtree.c
	$(CC) $(CFLAGS) -c spec.c
	$(CC) $(CFLAGS) butterfly.c -o butterfly quadtree.o spec.o

	$(CC) $(CFLAGS) bf_one_block.c -o bf_one_block quadtree.o
	$(CC) $(CFLAGS) build_quadtree.c -o build_quadtree quadtree.o
