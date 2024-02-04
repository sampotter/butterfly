Hierarchical Matrices
=====================

Concepts
--------

Transposition
^^^^^^^^^^^^^

First, some general principles about how transposition is handled in ``butterfly``:

- By default, "the transpose" refers to the Hermitian transpose. This is the same behavior as the ``'`` operator in MATLAB, but different from NumPy's ``.T`` property. That is, for a real matrix, the matrix is simply transposed; for a complex matrix, the matrix is transposed and its components are conjugated. Unfortunately, NumPy lacks something like a ``.H`` property (although this might change at some point---see `this issue <https://github.com/numpy/numpy/issues/13797>`_). Our view is that the Hermitian transpose is nearly always what is desired when working with complex matrices. In the rare case that the "true" transpose of a complex matrix is needed, and exception can be carved out.

Implementing consistent behavior for transposition in a hierarchical matrix library is tricky. Recall that the transpose of a block matrix satisfies:

.. math::
   \begin{bmatrix}
     A_{11} & \cdots & A_{1n} \\
	 \vdots & \ddots & \vdots \\
	 A_{m1} & \cdots & A_{mn}
   \end{bmatrix}^* = \begin{bmatrix}
     A_{11}^* & \cdots & A_{m1}^* \\
	 \vdots   & \ddots & \vdots   \\
	 A_{1n}^* & \cdots & A_{mn}^*
   \end{bmatrix}

The transpose of a hierarchical block matrix follows by recursive application of this equation. This is simple enough mathematically, but deciding what to do when a function like ``bfMatTranspose`` is applied to an argument with a hierarchical block structure is complicated.
