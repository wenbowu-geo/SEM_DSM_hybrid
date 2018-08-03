# -*- Autoconf -*-


## --------------------------- ##
## Autoconf macros for OPENCL. ##
## --------------------------- ##


# ----------------------------------------------------------------------
# CIT_OPENCL_CONFIG
# ----------------------------------------------------------------------
# Determine the directory containing <CL/cl.h>
AC_DEFUN([CIT_OPENCL_CONFIG], [
  AC_ARG_VAR(OCL_INC, [Location of OpenCL include files])
  AC_ARG_VAR(OCL_LIB, [Location of OpenCL library libOpenCL])

  AC_LANG_PUSH([C])
  AC_REQUIRE_CPP
  CFLAGS_save="$CFLAGS"
  LDFLAGS_save="$LDFLAGS"
  LIBS_save="$LIBS"

  dnl Check for OpenCL headers
  if test "x$OCL_INC" != "x"; then
    OCL_CFLAGS="-I$OCL_INC"
    CFLAGS="$OCL_CFLAGS $CFLAGS"
  fi
  # tests for Apple Mac OsX OpenCL header file
  case $host in
  *apple*) HEADER_H="OpenCL/cl.h";;
  *) HEADER_H="CL/cl.h";;
  esac
  AC_CHECK_HEADER([${HEADER_H}], [], [
    AC_MSG_ERROR([OpenCL header not found; try setting OCL_INC.])
  ])

  if test "x$OCL_LIB" != "x"; then
    # tests for Apple Mac OsX which could add "-framework OpenCL"
    case $host in
    *apple*) OCL_LDFLAGS="$OCL_LIB";;
    *) OCL_LDFLAGS="-L$OCL_LIB";;
    esac
    LDFLAGS="$OCL_LDFLAGS $LDFLAGS"
  else
    # adds default OpenCL library
    case $host in
    *apple*) OCL_LIBS="-framework OpenCL";;
    *) OCL_LIBS="-lOpenCL";;
    esac
  fi

  LIBS="$OCL_LIBS $LIBS"
  AC_MSG_CHECKING([for clGetPlatformIDs in OpenCL])
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
    ]],
  	             [[clGetPlatformIDs(0, 0, 0);]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([OpenCL library not found; try setting OCL_LIB.])
  ])

  CFLAGS="$CFLAGS_save"
  LDFLAGS="$LDFLAGS_save"
  LIBS="$LIBS_save"
  AC_LANG_POP([C])

  AC_SUBST([OCL_CFLAGS])
  AC_SUBST([OCL_LDFLAGS])
  AC_SUBST([OCL_LIBS])
])dnl CIT_OPENCL_COMPILER


dnl end of file
