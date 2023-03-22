# Guidelines

## Header files

1. **Put each function declaration on a single line, regardless of how long it is.** In Emacs compilation mode, if the declarations are split onto multiple lines, not all function arguments will appear in the compilation buffer.

## API

1. Most of the types in this library aren't opaque. This could be changed with a bit of complication, and likely at the cost of some efficiency. Leaving this point aside, this means we can directly malloc the types and assemble them ourselves. Nevertheless, it's preferable to define a `New` function for each type. This will allow us to use a different malloc at a later time easily, or to use more involved approaches to memory allocation. Also, if we do elect to make the types opaque later, then the change will be less painful.
