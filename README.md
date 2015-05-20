treatSens
=========

Utilities to investigate sensitivity to unmeasured confounding in parametric models with either binary or continuous treatment.

A package for R, with C/C++.

Pre-built binaries of the package are available on http://cran.r-project.org/web/packages/treatSens/index.html. These can be installed from within R using the typical `install.packages()` mechanism.

Steps to install from source:

  1. Install development tools for your operating system:

    1. Linux/Unix should already have this installed
    2. OS X:
        1. Xcode (https://developer.apple.com/xcode/downloads/)
        2. gfortran (https://gcc.gnu.org/wiki/GFortranBinaries#MacOS)
    3. Windows: Rtools (http://cran.r-project.org/bin/windows/Rtools/)

  2. Install the devtools package from within R:

    `install.packages("devtools")`

  3. Run:

    `install_github("vdorie/treatSens")`

Prior dependencies may have to be installed separately. A simple way to do this is to run the command:

    `install.packages(c("lme4", "dbarts"))`