CFLAGS := ${CFLAGS} -g -O0 -std=gnu99 -DBF_DEBUG -DBF_DOUBLE -lm
# CFLAGS := ${CFLAGS} -O3 -std=gnu99 -DBF_DOUBLE -lm -flto
CFLAGS := ${CFLAGS} -Wall -Wextra -Wshadow -pedantic
CFLAGS := ${CFLAGS} -lopenblas -I/usr/include/openblas

OBJECTS := dtype.o error.o fac.o geom.o helm2.o mat.o ptr_array.o	\
	quadtree.o rand.o splitmix64.o vec.o xoshiro256plus.o

all:
	$(CC) $(CFLAGS) -c dtype.c
	$(CC) $(CFLAGS) -c error.c
	$(CC) $(CFLAGS) -c fac.c
	$(CC) $(CFLAGS) -c geom.c
	$(CC) $(CFLAGS) -c helm2.c
	$(CC) $(CFLAGS) -c mat.c
	$(CC) $(CFLAGS) -c ptr_array.c
	$(CC) $(CFLAGS) -c quadtree.c
	$(CC) $(CFLAGS) -c rand.c
	$(CC) $(CFLAGS) -c splitmix64.c
	$(CC) $(CFLAGS) -c vec.c
	$(CC) $(CFLAGS) -c xoshiro256plus.c

	$(CC) $(CFLAGS) butterfly.c -o butterfly $(OBJECTS)
	$(CC) $(CFLAGS) bf_one_block.c -o bf_one_block $(OBJECTS)
