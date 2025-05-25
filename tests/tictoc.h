/*
 * tests/tictoc.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains function used for timing.
 */

#ifndef _TIMER_H_
#define _TIMER_H_

#include "config.h"

#ifdef _WIN32

#include <windows.h>
typedef DWORD tictoc_t;
#else

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
typedef struct timeval tictoc_t;
#else
#include <ctime>
typedef time_t tictoc_t;
#endif

#endif

void   tic(tictoc_t *tv);   /* start timing. */
double toc(tictoc_t *tv);   /* stop  timing. */

#endif  /* _TIMER_H_ */

