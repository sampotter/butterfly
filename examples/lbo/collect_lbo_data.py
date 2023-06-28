#!/usr/bin/env python

import argparse
import itertools as it
import numpy as np
import subprocess
import time

from pathlib import Path

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mesh_dir', type=str)

    args = parser.parse_args()

    mesh_paths = list(Path(args.mesh_dir).glob('*.obj'))



    results_dir_path = Path('results')
    results_dir_path.mkdir(exist_ok=True)

    results_raw_dir_path = results_dir_path/'raw'
    results_raw_dir_path.mkdir(exist_ok=True)
