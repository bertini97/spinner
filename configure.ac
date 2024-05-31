AC_INIT([spinner],[1.0],[lorenzobertini97@gmail.com])
AM_INIT_AUTOMAKE([gnu -Wall -Werror])

AC_PROG_CC
AM_PROG_AR
LT_INIT

AC_CONFIG_FILES([makefile])
AC_CHECK_HEADERS([stdint.h])
AC_CONFIG_MACRO_DIRS([m4])

AC_ARG_ENABLE([release],
    [AS_HELP_STRING([--enable-release], [build with optimization flags])],
    [enable_release=$enableval], [enable_release=no])

AM_CONDITIONAL([RELEASE_BUILD], [test "x$enable_release" = "xyes"])

AC_TYPE_SIZE_T
AC_FUNC_MALLOC

AC_CHECK_FUNCS([memset pow sqrt])

AC_OUTPUT
