#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT(treatSens, 3.0, vdorie@gmail.com)
AC_CONFIG_SRCDIR([src/R_interface.cpp])
AC_CONFIG_AUX_DIR([tools/build-aux])
AC_CONFIG_HEADER([src/config.hpp])
AC_CONFIG_HEADER([src/misc/config.h])
AC_CONFIG_HEADER([src/external/config.h])

m4_include([tools/m4/ax_check_compile_flag.m4])
m4_include([tools/m4/ax_pthread.m4])
m4_include([tools/m4/ax_compiler_vendor.m4])
m4_include([tools/m4/ax_compiler_ext.m4])

# Checks for programs.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CXXPICFLAGS=`"${R_HOME}/bin/R" CMD config CXXPICFLAGS`

AC_PROG_CC
AC_PROG_CC_C99
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP

AC_LANG(C)
AC_LANG_PUSH(C++)
AC_C_RESTRICT
AC_LANG_POP(C++)

AC_C_BIGENDIAN

AC_CANONICAL_HOST
AX_COMPILER_VENDOR
AC_SUBST(COMPILER_VENDOR, $ax_cv_c_compiler_vendor)
AX_COMPILER_EXT

# on SOLARIS, need to strip out some hardware flags from the library,
# so we test for those that are explicitly added with intrinsics
# but not supported by the CPU
mapfile_flag=
case $host_os in
  solaris*)
    case $host_cpu in
      i[[3456]]86*|x86_64*|amd64*)
        if test "$ax_cv_c_compiler_vendor" = sun; then
          m4_include([tools/m4/ax_ext_solaris.m4])
          m4_include([tools/m4/ax_gcc_x86_cpuid.m4])
          m4_include([tools/m4/ax_gcc_x86_avx_xgetbv.m4])
          AX_EXT_SOLARIS
          
          unsupported_isets=
          # only test for those that we explicitly add beyond what
          # the compiler will come up with on "-xarch=native"
          for iset in "avx;AVX" "sse4_1;SSE4_1" "avx2;AVX2"; do
            iset_varname=$(echo $iset | cut -d ";" -f 1)
            iset_mapname=$(echo $iset | cut -d ";" -f 2)
            # have variables test for CPU and OS support
            if test x"$(eval echo \$ax_cv_have_${iset_varname}_ext)" = x"no"; then
              unsupported_isets="${unsupported_isets} ${iset_mapname}"
            fi
          done
          if test ! -z "${unsupported_isets}"; then
            unsupported_isets=${unsupported_isets##[[:space:]]} # remove leading space
            echo "creating mapfile to remove ${unsupported_isets} instruction set(s)"
            echo "\$mapfile_version 2\nCAPABILITY {\n  HW -= ${unsupported_isets};\n};\n" > src/mapfile.solaris
            mapfile_flag="-M mapfile.solaris"
          fi
        fi
      ;;
    esac
  ;;
esac

AC_SUBST(MAPFILE_FLAG, $mapfile_flag)

# Checks for libraries.
AX_PTHREAD
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)

# Checks for header files.
AC_CHECK_HEADERS([sys/sysctl.h])
AC_LANG_PUSH(C++)
AC_CHECK_HEADERS([cstdint])
AC_LANG_POP(C++)

if test “$ac_cv_header_cstdint” == no
then
  AC_CHECK_HEADER([stdint.h])
fi

AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([malloc.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_CHECK_SIZEOF([size_t])
AC_CHECK_ALIGNOF([void*])

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([clock_gettime])
AC_CHECK_FUNCS([ffs])
AC_FUNC_ALLOCA

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
