# To-do list

1. [ ] `CamelCase` -> `snake_case`
2. [ ] more abbreviations (e.g., `BfMatDenseComplex` -> `bf_zmat`)
3. [ ] add SVDs to the matrix type hierarchy... ditto IDs eventually
4. [ ] add casts for handles (used in `*Delete`)
5. [ ] x-macros for automatically generating casts
6. [ ] functions in interfaces which aren't specialized should just be functions for the base type (e.g., `bfMatInstanceOf` shouldn't be in the interface... just do `bfMatInstanceOf(ToMat(...), matType)`
7. [ ] use `_Generic` to add generic casts, such as `to_mat`
8. [ ] should probably make quadtree uniform? see quadtree and points plots...
