CFLAGS := ${CFLAGS} -g -O0 -std=gnu99 -DBF_DEBUG -DBF_DOUBLE -lm
# CFLAGS := ${CFLAGS} -O3 -std=gnu99 -DBF_DOUBLE -lm -flto
CFLAGS := ${CFLAGS} -Wall -Wextra -Wshadow -pedantic
CFLAGS := ${CFLAGS} -lopenblas -I/usr/include/openblas

OBJECTS := dtype.o error.o fac.o geom.o helm2.o mat.o mat_block_coo.o	\
	mat_dense_complex.o mat_diag_real.o ptr_array.o quadtree.o rand.o	\
	splitmix64.o util.o vec.o xoshiro256plus.o

all:
	$(CC) $(CFLAGS) -c dtype.c
	$(CC) $(CFLAGS) -c error.c
	$(CC) $(CFLAGS) -c fac.c
	$(CC) $(CFLAGS) -c geom.c
	$(CC) $(CFLAGS) -c helm2.c
	$(CC) $(CFLAGS) -c mat.c
	$(CC) $(CFLAGS) -c mat_block_coo.c
	$(CC) $(CFLAGS) -c mat_block_dense.c
	$(CC) $(CFLAGS) -c mat_dense_complex.c
	$(CC) $(CFLAGS) -c mat_diag_real.c
	$(CC) $(CFLAGS) -c ptr_array.c
	$(CC) $(CFLAGS) -c quadtree.c
	$(CC) $(CFLAGS) -c rand.c
	$(CC) $(CFLAGS) -c splitmix64.c
	$(CC) $(CFLAGS) -c util.c
	$(CC) $(CFLAGS) -c vec.c
	$(CC) $(CFLAGS) -c xoshiro256plus.c

	$(CC) $(CFLAGS) butterfly.c -o butterfly $(OBJECTS)
	$(CC) $(CFLAGS) bf_one_block.c -o bf_one_block $(OBJECTS)
	$(CC) $(CFLAGS) bf_all_blocks.c -o bf_all_blocks $(OBJECTS)
