# butterfly

Brief explanation of the code:

- `bf_one_block.c` is a driver which tests out building a single
  butterfly factorization and writes the results to disk
  - run `make` to make it
  - run `python plot_bf_one_block.py` to make plots (will be saved to
    `./plots`)
- `fac.c` is where most of the action happens
- ignore `butterfly.c`
- `mat.c` contains a wrapper around OpenBLAS
- `helm2.c` has 2D Helmholtz-specific stuff
