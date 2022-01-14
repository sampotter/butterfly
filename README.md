# butterfly

Brief explanation of the code:

- run `make` to compile all the C code (this is a dumb and
  straightforward Makefile)
- `bf_one_block.c` is a driver which tests out building a single
  butterfly factorization and writes the results to disk
  - run `./make_test_data.py` followed by `./bf_one_block
    circle_with_16384_points.bin` to test it
  - run `python plot_bf_one_block.py` to make plots (will be saved to
    `./plots`)
- `bf_all_blocks.c` contains the start of a multilevel butterfly
  factorization (WIP!)
  - run `./bf_all_blocks circle_with_16384_points.bin` to test
  - run `./bf_plot_blocks.py` to plot hierarchical blocks and their
    levels
- `fac.c` is where most of the action happens
- ignore `butterfly.c`
- `helm2.c` has 2D Helmholtz-specific stuff

The files that start with `mat` contain the start of an
object-oriented matrix library. This is still in progress, but the
goal is to have a fairly generic algorithm which will take a
specification of a butterfly factorization and produce a hierarchical
block matrix implemented using `mat`. At the same time, this matrix
library should be flexibleenough that we can start experimenting with
direct solvers once the multilevel butterfly factorization works.
