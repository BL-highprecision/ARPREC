#ifndef PSLQ1_H
#define PSLQ1_H

#include "tictoc.h"
#include "matrix.h"

#undef inline

/* Mode codes */
#define NR_MODES                   1
#define MODE_ALGEBRAIC_TEST        0

/* Result codes */
#define RESULT_CONTINUE            0
#define RESULT_RELATION_FOUND      1
#define RESULT_PRECISION_EXHAUSTED 2

#define RESULT_RANGE_TOO_LARGE     3
#define RESULT_LARGE_VALUE         4
#define RESULT_VERY_LARGE_VALUE    5

#define RESULT_DUPLICATE           6

/* Constants */
#define LOG_2_BASE_10              3.01029995663981e-01
#define DEFAULT_GAMMA              1.154700538379252

/* Timer index constants */
#define NR_TIMERS                  7
#define TIMER_MP_UPDATE            0
#define TIMER_MPM_UPDATE           1
#define TIMER_MPM_LQ               2
#define TIMER_MP_INIT              3
#define TIMER_MPM_INIT             4
#define TIMER_PSLQ_TOTAL           5
#define TIMER_MPM_ITERATE          6

/* Timer macros */
#define TIMER_BEGIN(n)  { tictoc_t tv;  tic(&tv); 
#define TIMER_END(n)      timers(n) += toc(&tv); }

/* Precision control macros */
#define SET_PREC(n) mp::mpsetprecwords(n)
#define PREC_START  int old_nw = 0; \
                    if (nr_words) { \
                      old_nw = mp::mpgetprecwords(); \
                      mp::mpsetprecwords(nr_words); \
                    }
#define PREC_END    if (nr_words) { \
                      mp::mpsetprecwords(old_nw); \
                    }

/* Global variables */
extern int debug_level;
extern int pslq_iter;
extern double timers[];

/* Level 1 routines.  Found in pslq1_util.cpp and pslq1_templates.cpp. */
void clear_timers();
void report_timers();
void init_data(int mode, int n, int r, int s, 
               matrix<mp_real> &x, matrix<mp_real> &ans);

int reduce_pslq(matrix<mp_real> &h, matrix<mp_real> &y, 
                matrix<mp_real> &b, const mp_real &eps);
void init_pslq(const matrix<mp_real> &x, matrix<mp_real> &y, 
               matrix<mp_real> &h, matrix<mp_real> &b);
int iterate_pslq(double gamma, matrix<mp_real> &y, matrix<mp_real> &h, 
                 matrix<mp_real> &b, const mp_real &eps, const mp_real &teps);
int pslq1(const matrix<mp_real> &x, matrix<mp_real> &rel, const mp_real &eps, 
          double gamma = DEFAULT_GAMMA); 

/* Swaps the two elements x and y. */
/*
template <class T>
inline void swap(T &x, T &h) {
  T t = x;
  x = y;
  y = t;
}
*/

/* Swaps the two elementx x and y. 
   Specialization for mp_real type. */
/*
inline void swap(mp_real &x, mp_real &y) {
  mp_real::swap(x, y);
}
*/

template <class T>
int  matrix_minmax(const matrix<T> &v, T &v_min, T &v_max, bool isMin=true);
template <class T>
void lq_decomp(int n, int m, matrix<T> &a);
template <class T>
void matmul_left(const matrix<T> &a, matrix<mp_real> &b);
template <class T> 
void bound_pslq(const matrix<T> &h, T &r);

#include "pslq1_templates.cpp"

#endif

