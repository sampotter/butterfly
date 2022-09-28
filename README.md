# butterfly

*Overview to come...*

## Compilation

Use [Meson](https://mesonbuild.com/) to compile this library:
```
meson setup builddir
cd builddir
meson compile
```
This will build all of the examples, as well. Afterwards, the compiled executables for the examples will be in `./builddir/examples`.

### Error handling

This library features "OpenGL-style" error handling (e.g., [see this page](https://www.khronos.org/opengl/wiki/OpenGL_Error)). The guiding principles are three-fold:

1. It should allow you to pinpoint exactly where runtime user error occurs, useful for debugging the library.
2. If an error occurs when calling a function, the function should be a no-op; control should be procede normally afterwards (modulo the downstream effects of the error!).
3. It exists "in parallel" with the rest of the code, so that the error handling system can be completely stripped from a build if desired.

As this library is still in development, these principles are not fully realized yet, but this is the goal.

See [Intel Embree](https://www.embree.org/) for another example of a library with this style of error handling.

## Using Emacs lsp-mode

Emacs's [lsp-mode](https://emacs-lsp.github.io/lsp-mode/tutorials/CPP-guide/) can be used to add some IDE-like features to Emacs. For this to work with something like [clangd](https://clangd.llvm.org/) (recommended), clangd needs to be able to find a `compile_commands.json` file which describes the build. Meson generates this automatically, but it stores it in the build directory (e.g., `builddir` above). When opening a project file for the first time, lsp-mode will ask for the location of the "project root". Make sure that this directory contains a copy of the most recent `compile_commands.json` file.

## The hierarchical matrix types

This library includes a set of types for modeling recursively composed hierarchical matrices. These types support runtime polymorphism implemented using macros defined in [interface.h](./include/bf/interface.h). The root "class" in the hierarchcy is [BfMat](./include/bf/mat.h).
