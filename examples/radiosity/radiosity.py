#!/usr/bin/env python

import sys
sys.path.insert(-1, '..') # for util
sys.path.insert(-1, '../../wrappers/python') # for butterfly

import butterfly as bf

if __name__ == '__main__':
    trimesh = bf.Trimesh.from_obj('67p.obj')
    F = bf.MatCsrReal.new_view_factor_matrix_from_trimesh(trimesh)
    print(F)
