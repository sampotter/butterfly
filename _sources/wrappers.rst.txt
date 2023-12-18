Wrappers
========

The goal of this library is to provide wrappers for the high-level languages commonly used in computational science and engineering. Eventually, we aim to provide wrappers for (at least) Python, MATLAB, and Julia. Currently, only the Python bindings are under development.

The guiding principle for the wrappers is to try to map the object-oriented approach of ``butterfly`` as naturally as possible onto an low-level library implemented in the "native style" of the high-level language. You can contrast this with the approach taken by gmsh, where an API generator is used to spit out bindings for as many high-level languages as possible. Our view is that, while API generators are seductive, hand-written bindings which more naturally map onto the idiom of the high-level language are likely to be more useful in the long run. In particular, the goal should be to have the low-level wrappers be useful enough to work with directly, rather than serve as a substrate for a higher-level library.

Python
------

Several points:

- The types in ``butterfly`` map directly onto Python extension types, defined using Cython's ``cdef class``.

Python's weak duck typing means there is no need to cast between different types like in strongly typed languages like C++ and, indeed, in the object system used by ``butterfly``. There's a tension here: since ``butterfly`` types map onto Python extension types, we need to make sure that we always return the most derived type from C to Python. To this end, we refer to the process of casting a C ``butterfly`` type to its "real" (most derived) type and instantiating and returning the corresponding extension type as "reification".

For example, ``reify_mat`` takes a ``BfMat *``. If the real type of the argument is ``BfMatDenseComplex``, we use ``bfMatToMatDenseComplex`` to cast to this type, wrap the result in ``MatDenseComplex``, and return that. Filling a container with a heterogenous collection of instances which all share a common base type is no problem.

MATLAB
------

Julia
-----
