/* include/arprec/arprec_config.h.  Generated from arprec_config.h.in by configure.  */
#ifndef ARPREC_CONFIG
#define ARPREC_CONFIG
/*
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001-02
 *
 * Contains configuration options used by the ARPREC headers. 
 */

/* Define to 1 to enable debugging code. */
#ifndef ARPREC_DEBUG
/* #undef ARPREC_DEBUG */
#endif

/* Define to 1 to use valgrind debugging code. Requires ARPREC_DEBUG to be
   set. */
#ifndef ARPREC_VALGRIND
/* #undef ARPREC_VALGRIND */
#endif

/* Define to 1 to inline commonly used functions. */
#ifndef ARPREC_INLINE
#define ARPREC_INLINE 1
#endif

/* Define to something appropriate if arprec APIs needs special
   flag in order to be exported (needed in MSVC when making DLL).  */
#ifndef ARPREC_API
#define ARPREC_API 
#endif

/* Define to 1 to use the quad-double library. */
#ifndef ARPREC_QD
/* #undef ARPREC_QD */
#endif

#endif
