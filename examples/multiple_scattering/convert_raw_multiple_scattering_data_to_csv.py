#!/usr/bin/env python

import csv

from glob import glob
from pathlib import Path

def get_row_from_raw_data_file(path):
    k, r, a, b, h = [float(_[1:]) for _ in str(path).split('/')[-1][:-4].split('_')]
    with open(path, 'r') as f:
        lines = f.readlines()

        line = next(_ for _ in lines if 'points per wavelength' in _)
        ppw = float(line.split()[-1])

        line = next(_ for _ in lines if 'ellipses [' in _)
        ne = int(line.split()[-3])

        line = next(_ for _ in lines if 'points [' in _)
        nx = int(line.split()[-3])

        line = next(_ for _ in lines if 'Assembling butterflied' in _)
        bf_T_assemble = float(line.split()[-1][1:-2])

        line = next(_ for _ in lines if 'Solving butterflied system using preconditioned GMRES' in _)
        bf_pGMRES_nit = int(line.split()[-3])
        bf_pGMRES_T = float(line.split()[-1][1:-2])

        line = next(_ for _ in lines if 'Solving FMM system using preconditioned GMRES' in _)
        fmm_pGMRES_nit = int(line.split()[-3])
        fmm_pGMRES_T = float(line.split()[-1][1:-2])

        line = next(_ for _ in lines if 'relative error between butterfly and FMM MVPs' in _)
        MVP_rel_err_bf_vs_FMM = float(line.split()[-1])

        line = next(_ for _ in lines if 'rel l2 error in sigma (FMM (preconditioned GMRES) vs butterfly (preconditioned GMRES))' in _)
        sigma_rel_err_bf_vs_FMM = float(line.split()[-1])

    return {
        'k': k,
        'r': r,
        'a': a,
        'b': b,
        'h': h,
        'ppw': ppw,
        'ne': ne,
        'nx': nx,
        't (BF assemb)': bf_T_assemble,
        'nit (BF pGMRES)': bf_pGMRES_nit,
        't (BF pGMRES)': bf_pGMRES_T,
        'nit (FMM pGMRES)': fmm_pGMRES_nit,
        't (FMM pGMRES)': fmm_pGMRES_T,
        'MVP rel err (bf vs FMM)': MVP_rel_err_bf_vs_FMM,
        'sigma rel err (bf vs FMM)': sigma_rel_err_bf_vs_FMM}

if __name__ == '__main__':
    results_dir_path = Path('results')
    assert results_dir_path.exists()

    results_raw_dir_path = results_dir_path/'raw'
    assert results_raw_dir_path.exists()

    raw_data_file_paths = results_raw_dir_path.glob('*.txt')

    csv_path = results_dir_path/'multiple_scattering_data.csv'
    assert not csv_path.exists()

    fieldnames = ['k', 'r', 'a', 'b', 'h', 'ppw', 'ne', 'nx', 't (BF assemb)', 'nit (BF pGMRES)', 't (BF pGMRES)', 'nit (FMM pGMRES)', 't (FMM pGMRES)', 'MVP rel err (bf vs FMM)', 'sigma rel err (bf vs FMM)']

    with open(csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        csv_writer.writeheader()
        for raw_data_file_path in raw_data_file_paths:
            print(raw_data_file_path)
            row = get_row_from_raw_data_file(raw_data_file_path)
            csv_writer.writerow(row)
