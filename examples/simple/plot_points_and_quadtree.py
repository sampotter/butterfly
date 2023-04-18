#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

def read_quadtree_boxes(path):
    with open(path, 'r') as f:
        for line in f.readlines():
            yield map(float, line.split())

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'usage: {sys.argv[0]} QUADTREE.TXT POINTS.BIN')
        exit(1)

    plt.figure()

    quadtree_boxes = []
    for x0, x1, y0, y1 in read_quadtree_boxes(sys.argv[1]):
        rect = Rectangle((x0, y0), x1 - x0, y1 - y0)
        quadtree_boxes.append(rect)
    pc = PatchCollection(quadtree_boxes, facecolor='none', edgecolor='k')
    plt.gca().add_collection(pc)

    x, y = np.fromfile(sys.argv[2]).reshape(-1, 2).T

    plt.scatter(x, y, s=1, c='r', zorder=2)

    plt.gca().set_aspect('equal')

    plt.show()
