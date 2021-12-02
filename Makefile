CFLAGS := ${CFLAGS} -g -O0 -std=gnu99 -DBF_DEBUG -DBF_DOUBLE -lm -Wall -Wextra
CFLAGS := ${CFLAGS} -lopenblas -I/usr/include/openblas

all:
	$(CC) $(CFLAGS) -c ptr_array.c
	$(CC) $(CFLAGS) -c dtype.c
	$(CC) $(CFLAGS) -c geom.c
	$(CC) $(CFLAGS) -c helm2.c
	$(CC) $(CFLAGS) -c mat.c
	$(CC) $(CFLAGS) -c quadtree.c
	$(CC) $(CFLAGS) -c rand.c
	$(CC) $(CFLAGS) -c spec.c
	$(CC) $(CFLAGS) -c splitmix64.c
	$(CC) $(CFLAGS) -c xoshiro256plus.c

	$(CC) $(CFLAGS) butterfly.c -o butterfly ptr_array.o dtype.o geom.o \
		helm2.o mat.o quadtree.o rand.o spec.o splitmix64.o xoshiro256plus.o
	$(CC) $(CFLAGS) bf_one_block.c -o bf_one_block ptr_array.o dtype.o geom.o \
		helm2.o mat.o quadtree.o rand.o spec.o splitmix64.o xoshiro256plus.o
	$(CC) $(CFLAGS) build_quadtree.c -o build_quadtree ptr_array.o dtype.o \
		geom.o helm2.o mat.o rand.o quadtree.o spec.o splitmix64.o \
		xoshiro256plus.o
