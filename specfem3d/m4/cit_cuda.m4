# -*- Autoconf -*-


## ------------------------- ##
## Autoconf macros for CUDA. ##
## ------------------------- ##


# ----------------------------------------------------------------------
# CIT_CUDA_CONFIG
# ----------------------------------------------------------------------
# Determine the directory containing <cuda_runtime.h>
AC_DEFUN([CIT_CUDA_CONFIG], [

  # influential environment variables
  AC_ARG_VAR(NVCC, [NVIDIA CUDA compiler command])
  AC_ARG_VAR(CUDA_FLAGS, [CUDA compiler flags])
  AC_ARG_VAR(CUDA_INC, [Location of CUDA include files])
  AC_ARG_VAR(CUDA_LIB, [Location of CUDA library libcudart])

  # tests NVCC variable
  AS_IF([test x"$NVCC" = x],[
    NVCC=nvcc
  ])

  # Check for compiler
  # checks if program in path
  AC_PATH_PROG(NVCC_PROG, $NVCC)
  if test -z "$NVCC_PROG" ; then
    AC_MSG_ERROR([cannot find '$NVCC' program, please check your PATH.])
  fi

  # Checks for compiling and linking
  AC_LANG_PUSH([C])
  AC_REQUIRE_CPP
  CFLAGS_save="$CFLAGS"
  LDFLAGS_save="$LDFLAGS"
  LIBS_save="$LIBS"

  # uses nvcc compiler
  CFLAGS="$CUDA_FLAGS"
  if test "x$CUDA_INC" != "x"; then
    CUDA_CPPFLAGS="-I$CUDA_INC"
    CFLAGS="$CFLAGS $CUDA_CPPFLAGS"
  fi

  # Check for CUDA headers
  # runs test with nvcc
  AC_MSG_CHECKING([for cuda_runtime.h])
  ac_compile='$NVCC -c $CFLAGS conftest.$ac_ext >&5'
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[
#include <cuda.h>
#include <cuda_runtime.h>]],[[void* ptr = 0;]])
  ], [
    AC_MSG_RESULT(yes)
  ], [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([CUDA runtime header not found; try setting CUDA_INC.])
  ])

  # Check fo CUDA library
  if test "x$CUDA_LIB" != "x"; then
    CUDA_LDFLAGS="-L$CUDA_LIB"
    LDFLAGS="$CUDA_LDFLAGS $LDFLAGS"
  fi
  CUDA_LIBS="-lcudart"
  LIBS="$CUDA_LIBS $LIBS"

  # runs compilation test with nvcc
  AC_MSG_CHECKING([nvcc compilation with cudaMalloc in -lcudart])
  ac_compile='$NVCC -c $CFLAGS conftest.$ac_ext >&5'
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[
#include <cuda.h>
#include <cuda_runtime.h>]],[[void* ptr = 0;cudaMalloc(&ptr, 1);]])
  ], [
    AC_MSG_RESULT(yes)
  ], [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([CUDA library function with nvcc compilation failed; try setting CUDA_INC.])
  ])

  # runs linking test with nvcc
  AC_MSG_CHECKING([nvcc linking with cudaMalloc in -lcudart])
  ac_link='$NVCC -o conftest$ac_exeext $CFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>]],[[void* ptr = 0;cudaMalloc(&ptr, 1);]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([CUDA library linking with nvcc failed; try setting CUDA_LIB.])
  ])

  # runs linking test with standard compiler
  AC_MSG_CHECKING([linking with cudaMalloc in -lcudart])
  # C compiler linking
  ac_link='$NVCC -c $CFLAGS conftest.$ac_ext >&5; $CC -o conftest$ac_exeext $LDFLAGS conftest.$ac_objext $LIBS >&5'
  # Fortran compiler linking
  # note: fortran linking with ifort would need -nofor-main to succeed
  #ac_link='$NVCC -c $CFLAGS conftest.$ac_ext >&5; $FC -o conftest$ac_exeext $LDFLAGS conftest.$ac_objext $LIBS >&5'
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>]],[[void* ptr = 0;cudaMalloc(&ptr, 1);]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([CUDA library linking failed; try setting CUDA_LIB.])
  ])

  CFLAGS="$CFLAGS_save"
  LDFLAGS="$LDFLAGS_save"
  LIBS="$LIBS_save"
  AC_LANG_POP([C])

  AC_SUBST([NVCC])
  AC_SUBST([CUDA_CPPFLAGS])
  AC_SUBST([CUDA_LDFLAGS])
  AC_SUBST([CUDA_LIBS])
])dnl CIT_CUDA_COMPILER


dnl end of file
