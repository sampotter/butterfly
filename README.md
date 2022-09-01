# butterfly

*Overview to come...*

## Compilation

Use Meson to compile this library:
```
meson setup builddir
cd builddir
meson compile
```
This will build all of the examples, as well. Afterwards, the compiled executables for the examples will be in `./builddir/examples`.

### Using Emacs lsp-mode

Emacs's [lsp-mode](https://emacs-lsp.github.io/lsp-mode/tutorials/CPP-guide/) can be used to add some IDE-like features to Emacs. For this to work with something like [clangd](https://clangd.llvm.org/) (recommended), clangd needs to be able to find a `compile_commands.json` file which describes the build. Meson generates this automatically, but it stores it in the build directory (e.g., `builddir` above). When opening a project file for the first time, lsp-mode will ask for the location of the "project root". Make sure that this directory contains a copy of the most recent `compile_commands.json` file.
