#ifndef PSLQ3_H
#define PSLQ3_H

#include "pslq1.h"
#include "pslq2.h"
#include "matrix.h"

int pslq3(const matrix<mp_real> &x, matrix<mp_real> &rel, 
          const mp_real &eps, double gamma = DEFAULT_GAMMA);

template <class T>
int check_size_pslq(const matrix<T> &a, const matrix<T> &b, 
                    const T &mx1, const T &mx2);

#include "pslq3_templates.cpp"

#endif
