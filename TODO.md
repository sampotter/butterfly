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

For reference counting, retain functions may not actually be necessary? Not sure.

# Don't get rid of const?

Re: "Get rid of const...? and "Handle ownership correctly", if we *don't* get rid of const, then one issue is views. E.g., if we get a view of some subblock of a matrix which we have a constant pointer to, that view should also be const. Realistically, there are two levels of const here: shallow const (`BfMat const *` vs `BfMat *`) and deep const (`BfMatViewConst` vs `BfMatView`---these don't exist at the time of writing...). I don't currently have the time to build out the type hierarchy to model the latter accurately. But at least for shallow const, if we get a view of a `BfMat const *` subblock, we'd like it to also be `BfMat const *`. The problem with this is that we need to allocate memory to get this view. Critically, `free` *does not* take a const pointer. One benefit of using our own wrappers for malloc and free is that we can sidestep this problem. But even better would be to outlaw `free` entirely by only using `retain` and `release` type functions (see "Handle ownership correctly"). In this case, it makes perfectly good sense to have `retain` and `release` take const pointers. They don't *logically* modify the value being pointed to, they only increment or decrement the associated reference count.

# Add a top type

# All objects should be passed as pointers

Since we rely so heavily on using `NULL` as a tombstone value, we really need to be consistent about passing around objects "by reference" only.

# Add a NodeSpan type

There's a lot of useful stuff in [fac_streamer.c](./src/fac_streamer.c) related to working with contiguous spans of nodes. Given how much stuff we have to do with node spans when working with butterfly factorizations, it seems useful to break this off into a separate ADT.

# Provenance tracking

If I write my own interface to malloc, I can wrap it with a macro which will let me track where every object was allocated (i.e., which line number, what CPU time... could annotate with various bits of data).

# Invalidating objects

Objects should always be put in an obvious invalid state immediately after being created, with all values set by default to "bad value" for the type.

Once this is set up, it should be possible to set a flag to strip all of this out when running in release mode.

# Better interface for constructing objects

Add `Alloc` function which does what `New` does now. Create `New*` versions of `Init*` functions which do the same thing but also allocate. Then `Alloc`:`Init`:`New`::`Dealloc`:`Deinit`:`Delete` works as an anlogy, which should make things a little more regular and intuitive.

# Less error-handling noise

Realistically, it isn't possible to recover from an OOM. We should wrap malloc and catch this there, if anything. Getting rid of `BF_ERROR_MEMORY_ERROR` would remove a significant amount of error handling line noise, leaving errors which are more significant.
