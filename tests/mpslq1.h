#ifndef MPSLQ1_H
#define MPSLQ1_H

#include <arprec/mp_real.h>
#include "pslq1.h"

int iterate_mpslq(double gamma, matrix<mp_real> &y, matrix<mp_real> &h, 
                  matrix<mp_real> &b, const mp_real &eps, const mp_real &teps, 
                  bool use_only_one_pair);
int mpslq1(const matrix<mp_real> &x, matrix<mp_real> &rel, 
           const mp_real &eps, double gamma = DEFAULT_GAMMA); 

template <class T>
void sort_pslq(const matrix<T> &x, matrix<int> &ip);

template <class T>
int check_history(const matrix<T> &y, const matrix<T> &hist, const T &teps);

template <class T> 
void insert_history(const matrix<T> &y, matrix<T> &hist);

#include "mpslq1_templates.cpp"

#endif
