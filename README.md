treatSens
=========

Utilities to investigate sensitivity to unmeasured confounding in parametric models with either binary or continuous treatment.

A package for R, with C/C++.

The pre-built binary on CRAN has been archived. Please install from source:

1. Install development tools for your operating system:
    1. Linux/Unix should already have this installed; if not, use your package manager to install a C/C++ compiler.
    2. OS X: [XCode and gfortran](https://mac.r-project.org/tools/)
    3. Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

2. Install the `remotes` package from within R:

```R
install.packages("remotes")
```

3. Run:

```R
remotes::install_github("vdorie/treatSens")
```
