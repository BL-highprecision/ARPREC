/*
 * include/arprec/small_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef _MP_SMALL_INLINE_H_
#define _MP_SMALL_INLINE_H_
/* This file contains short inline mathematical routines.
 *   
 *  Many of these routines rely on IEEE floating points
 *  numbers to preform mathematical 'tricks'.
 */

double quick_two_sum(double a, double b, double &err);
double quick_two_diff(double a, double b, double &err);
double two_sum(double a, double b, double &err);
double two_diff(double a, double b, double &err);
double two_sqr(double a, double &err);

double SLOPPY_ANINT_POSITIVE(double a);
double FLOOR_POSITIVE(double a);
double CEIL_POSITIVE(double a);
double AINT(double a);
double POSITIVE_AINT(double a);
double ANINT(double a);
double SLOPPY_ANINT(double a);

void dd_add_d(double a[], double b, double c[]);
void dd_add_dd(double a[], double b[], double c[]);

double mp_two_prod(double a, double b, double &err);
double mp_two_prod_positive(double a, double b, double &err);

#ifndef ARPREC_FMS
void split(double a, double &hi, double &lo);
#endif

#ifdef ARPREC_INLINE
#include "small_inline.cpp"
#endif

#endif
