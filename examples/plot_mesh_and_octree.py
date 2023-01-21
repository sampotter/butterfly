import argparse
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

mesh_path = 'dragon0.obj'
octree_boxes_path = 'octree_boxes.txt'

octree_boxes = dict()
with open(octree_boxes_path, 'r') as f:
    for line in f:
        tok = line.split()
        depth = int(tok[0])
        if depth not in octree_boxes:
            octree_boxes[depth] = list()
        octree_boxes[depth].append([float(_) for _ in tok[1:]])
print(f'parsed octree boxes')

max_depth = max(octree_boxes.keys())
print(f'{max_depth = }')

mesh = pv.read(mesh_path)

plotter = pvqt.BackgroundPlotter()
plotter.background_color = 'white'
plotter.add_mesh(mesh)
for depth in [5]:
    if depth in octree_boxes:
        for bounds in octree_boxes[depth]:
            box = pv.Box(bounds)
            edges = box.extract_all_edges()
            plotter.add_mesh(edges, color='black')
