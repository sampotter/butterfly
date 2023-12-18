Todo
====

1. ARPACK isn't reentrant. This poses two problems:

   a. Obviously, it will be a bottleneck when it comes time to parallelize butterfly.

   b. Even for serial code, there may be cases where we want to have several instances of ARPACK running simultaneously so that we can batch together MVPs (this will happen with deterministic fast direct solvers).
