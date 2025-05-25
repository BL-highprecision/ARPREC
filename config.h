/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 to use FPU fix whenever C interface is called. */
/* #undef ALWAYS_FIX_FPU */

/* Any special symbols needed for exporting APIs. */
#define ARPREC_API 

/* Define to 1 to enable debugging code. */
/* #undef ARPREC_DEBUG */

/* If fused multiply-add is available, define correct macro for using it. */
/* #undef ARPREC_FMA */

/* If fused multiply-subtract is available, define correct macro for using it.
   */
/* #undef ARPREC_FMS */

/* Define to 1 to inline commonly used functions. */
#define ARPREC_INLINE 1

/* Define to 1 to use the QD library. */
/* #undef ARPREC_QD */

/* Set to 1 if using VisualAge C++ compiler for __fmadd builtin. */
/* #undef ARPREC_VACPP_BUILTINS_H */

/* Define to 1 to use valgrind debugging code. Requires ARPREC_DEBUG to be
   set. */
/* #undef ARPREC_VALGRIND */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to alternate name for `main' routine that is called from a `main' in
   the Fortran libraries. */
#define FC_MAIN main

/* Define to 1 if you have the `copysign' function. */
#define HAVE_COPYSIGN 1

/* Define to 1 if you have the `drand48' function. */
#define HAVE_DRAND48 1

/* Define to 1 if Fortran interface is to be compiled. */
#define HAVE_FORTRAN 1

/* Define to 1 if you have the <fpu_control.h> header file. */
/* #undef HAVE_FPU_CONTROL_H */

/* Define to 1 if you have the <fp.h> header file. */
/* #undef HAVE_FP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* define if the compiler supports basic templates */
#define HAVE_TEMPLATES 

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* arprec major version number */
#define MAJOR_VERSION 2

/* arprec minor version number */
#define MINOR_VERSION 8

/* Name of package */
#define PACKAGE "arprec"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "dhbailey@lbl.gov"

/* Define to the full name of this package. */
#define PACKAGE_NAME "arprec"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "arprec 2.2.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "arprec"

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.2.1"

/* arprec patch number (sub minor version) */
#define PATCH_VERSION 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "2.2.1"

/* Whether to use x86 fpu fix. */
#define X86 1
