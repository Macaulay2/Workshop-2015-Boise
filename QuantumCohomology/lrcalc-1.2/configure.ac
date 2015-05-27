# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([lrcalc],[1.2],[asbuch at math rutgers edu])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu])
AM_MAINTAINER_MODE

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

LT_INIT(win32-dll)

dnl Copied from MPFR 3.1.1 and modified.
dnl Configs for Windows DLLs.
dnl libtool requires "-no-undefined" for win32 dll
dnl
dnl "-Wl,output-def" is used to get a .def file for use by MS lib to make
dnl a .lib import library, described in the manual.
AC_SUBST(LIBLRCALC_LDFLAGS)
case $host in
  *-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
    if test "$enable_shared" = yes; then
      LIBLRCALC_LDFLAGS="$LIBLRCALC_LDFLAGS -no-undefined -Wl,--output-def,.libs/liblrcalc-0.def"
    fi
  ;;
esac

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset strtol])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
