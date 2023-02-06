# To-do list

1. [ ] `CamelCase` -> `snake_case`
2. [ ] more abbreviations (e.g., `BfMatDenseComplex` -> `bf_zmat`)
3. [ ] add SVDs to the matrix type hierarchy... ditto IDs eventually
4. [ ] add casts for handles (used in `*Delete`)
5. [ ] x-macros for automatically generating casts
6. [ ] functions in interfaces which aren't specialized should just be functions for the base type (e.g., `bfMatInstanceOf` shouldn't be in the interface... just do `bfMatInstanceOf(ToMat(...), matType)`
7. [ ] use `_Generic` to add generic casts, such as `to_mat`
8. [ ] should probably make quadtree uniform? see quadtree and points plots...
9. [ ] don't use my own typedefs for `double` and `double _Complex`---simpler to reason about and paves the way for supporting single precision at the same time
10. [ ] need to come up with a smarter implementation of multiple dispatch

# New type names

- `BfMat` -> `bf_mat`
- `BfMatDenseComplex` -> `bf_zmat` (+ `bf_cmat` for single-precision)
- `BfMatDenseReal` -> `bf_dmat` (+ `bf_smat` for single-precision)
- `BfMatBlock` -> `bf_bmat` ??? (eliminate???)
- `BfMatBlockDense` -> `bf_bmat` ???
- `BfMatBlockCoo` -> `bf_bcoo`
- `BfMatBlockDiag` -> `bf_bdiag`
- `BfMatDiagReal` -> `bf_ddiag` (+ `bf_sdiag` for single-precision)
- `BfMatDiagComplex` -> `bf_zdiag` (+ `bf_cdiag` for single-precision)
- `BfMatProduct` -> `bf_matprod`
- `BfVec` -> `bf_vec`
- `BfVecReal` -> `bf_dvec` (+ `bf_svec` for single-precision)
- `BfVecComplex` -> `bf_zvec` (+ `bf_cvec` for single-precision)

and some new types:

- `bf_*coo`, w/ `*` in `{s, d, c, z}` for sparse COO matrix, e.g.
- `bf_*csr`, w/ `*` in `{s, d, c, z}` for sparse CSR matrix, e.g.

# Not implementing parts of the interfaces

Instead of using stubs (`BF_STUB(...)`), we should populate the
vtables with `NULL` function pointers by default. If a delegating
function encounters a `NULL` function pointer, it should just set the
error code and bail. This will reduce the amount of boilerplate we
need to write, but it will also make the intent of the code
clearer. It will be okay for functions to simply not override parts of
the interface. Especially since it's possible to dynamically update
the vtable in C (as compared to C++), this is fine. There are places
where it is better to not implement functions in the interface. For
instance, if we want to scale the rows or columns of a matrix, we can
actually modify the contents of the matrix, or we can pre- or
post-multiply by a diagonal matrix. Different behaviors may be
preferred under different circumstances.

# Debug vtable wrappers and instrumentation

For debugging or profiling, we may want to "wrap" a vtable. It
shouldn't be too hard to wrap particular functions in an interface in
order to intercept function calls. For instance, we might want to log
every memory allocation, or every matrix multiplication.

As another example, when it comes time to generate instructions for
our "matrix algebra VM" ("Maeve"), we may just want to replace matrix
vector product calls with functions which generate instructions and
record them somewhere.

# ~~Get rid of macros for interfaces~~ (Done!)

These aren't actually helpful. Using them, we tie ourself in a knot
which makes it hard to accomplish the sort of things described above
in "Debug vtable wrappers and instrumentation".

# Get rid of const...?

It's not clear how useful `const` actually is. If I try to be too fussy about `const` correctness in C, then it seems like the main downstream effect is having to maintain a parallel set of data structures, one for `const` pointers, and one more regular pointers.

The bigger problem is that `const` in C means "physical" `const`ness: if a variable is `const`, then we can't modify its bitwise representation. For example, if we have a const struct with a non-const pointer, we can still modify what is pointed to by that pointer. So, we do not get logical `const`ness from physical `const`ness.

What we want is logical `const`ness to aid with debugging and program correctness, but it does not seem possible to achieve this with `const` alone---it isn't a powerful enough tool for the job. Instead, it might be better to juts abandon the project and reap benefits elsewhere---for instance, there being less code.

The question is where or not it's helpful to use `const` anywhere beyond this...

# Useful libraries

- https://librsb.sourceforge.net/
- https://github.com/opencollab/arpack-ng

# Handle ownership correctly

As we build up more elaborate data types, we're going to need to think more carefully about how we handle the ownership of our different types.

One straightforward approach would be to use reference counting with explicit `retain` and `release` functions, a la Embree. Another option would be to have ownership semantics which are set manually.

# Add a top type

# All objects should be passed as pointers

Since we rely so heavily on using `NULL` as a tombstone value, we really need

# Add a NodeSpan type

There's a lot of useful stuff in [fac_streamer.c](./src/fac_streamer.c) related to working with contiguous spans of nodes. Given how much stuff we have to do with node spans when working with butterfly factorizations, it seems useful to break this off into a separate ADT.
