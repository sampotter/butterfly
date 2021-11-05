#pragma once

#include <stdlib.h>

typedef double point2[2];

typedef struct quadtree quadtree_s;
typedef struct quadtreeNode quadtreeNode_s;

struct quadtreeNode {
	quadtree_s const *root;
	size_t offset[5]; // offset[0] == 0 and offset[4] == size are sentinel values
	quadtreeNode_s *child[4];
};

struct quadtree {
	point2 const *points;
	size_t num_points;
	size_t *perm;
	quadtreeNode_s *root;
};
