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
