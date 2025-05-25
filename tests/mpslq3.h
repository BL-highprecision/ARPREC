#ifndef MPSLQ3_H
#define MPSLQ3_H

#include "mpslq2.h"
#include "pslq3.h"

template <class T>
int update_mpslq2(const matrix<T> &da, const matrix<T> &db, 
		              matrix<mp_real> &a, matrix<mp_real> &b, 
                  matrix<mp_real> &h, matrix<mp_real> &y, 
									const mp_real &eps, const mp_real &teps);

int mpslq3(const matrix<mp_real> &x, matrix<mp_real> &rel, 
           const mp_real &eps, double gamma = DEFAULT_GAMMA);

#include "mpslq3_templates.cpp"

#endif
