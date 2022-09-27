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
