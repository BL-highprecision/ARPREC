/*
 * include/arprec/c_mp.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 * Contains C wrapper function prototypes for multi-precision arithmetic.
 * This can also be used from fortran code.
 */
#ifndef ARPREC_C_MP_H
#define ARPREC_C_MP_H

#include <arprec/arprec_config.h>

#ifdef __cplusplus
extern "C" {
#endif

  /* add */
  ARPREC_API void c_mpadd(const double *a, const double *b, double *c);
  ARPREC_API void c_mpadd_d(const double *a, double b, double *c);
  ARPREC_API void c_mpadd_ji(const double *a, int b, double *c);
  ARPREC_API void c_mpadd_jd(const double *a, double b, double *c);
  ARPREC_API void c_mpadd_zq(const double *a, const double *b, double *c);
  ARPREC_API void c_mpadd_zx(const double *a, double br, double bi, double *c);
  ARPREC_API void c_mpadd_zz(const double *a, const double *b, double *c);

  /* sub */
  ARPREC_API void c_mpsub(const double *a, const double *b, double *c);
  ARPREC_API void c_mpsub_d(const double *a, double b, double *c);
  ARPREC_API void c_mpsub_dq(double a, const double *b, double *c);
  ARPREC_API void c_mpsub_ji(const double *a, int b, double *c);
  ARPREC_API void c_mpsub_ij(int a, const double *b, double *c);
  ARPREC_API void c_mpsub_jd(const double *a, double b, double *c);
  ARPREC_API void c_mpsub_dj(double a, const double *b, double *c);
  ARPREC_API void c_mpsub_zq(const double *a, const double *b, double *c);
  ARPREC_API void c_mpsub_qz(const double *a, const double *b, double *c);
  ARPREC_API void c_mpsub_zx(const double *a, double br, double bi, double *c);
  ARPREC_API void c_mpsub_xz(double ar, double ai, const double *b, double *c);
  ARPREC_API void c_mpsub_zz(const double *a, const double *b, double *c);
  ARPREC_API void c_mpneg_q(const double *a, double *c);
  ARPREC_API void c_mpneg_z(const double *a, double *c);

  /* mul */
  ARPREC_API void c_mpmul(const double *a, const double *b, double *c);
  ARPREC_API void c_mpmul_ji(const double *a, int b, double *c);
  ARPREC_API void c_mpmul_qd(const double *a, double b, double *c);
  ARPREC_API void c_mpmul_qi(const double *a, int b, double *c);
  ARPREC_API void c_mpmul_zq(const double *a, const double *b, double *c);
  ARPREC_API void c_mpmul_zz(const double *a, const double *b, double *c);
  ARPREC_API void c_mpmul_zd(const double *a, double b, double *c);

  /* div */
  ARPREC_API void c_mpdiv(const double *a, const double *b, double *c);
  ARPREC_API void c_mpdiv_jj(const double *a, const double *b, double *c);
  ARPREC_API void c_mpdiv_ji(const double *a, int b, double *c);
  ARPREC_API void c_mpdiv_ij(int a, const double *b, double *c);
  ARPREC_API void c_mpdiv_qi(const double *a, int b, double *c);
  ARPREC_API void c_mpdiv_iq(int a, const double *b, double *c);
  ARPREC_API void c_mpdiv_qd(const double *a, double b, double *c);
  ARPREC_API void c_mpdiv_dq(double a, const double *b, double *c);
  ARPREC_API void c_mpdiv_zq(const double *a, const double *b, double *c);
  ARPREC_API void c_mpdiv_qz(const double *a, const double *b, double *c);
  ARPREC_API void c_mpdiv_zz(const double *a, const double *b, double *c);
  ARPREC_API void c_mpdiv_zd(const double *a, double b, double *c);
  ARPREC_API void c_mpdiv_dz(double a, const double *b, double *c);
  ARPREC_API void c_mpdmc(double a, double *b);
  ARPREC_API void c_mpmdc(const double *a, double *b, int *n);

  /*Special Functions*/
  ARPREC_API void c_mperfc(const double *a, double *b);
  ARPREC_API void c_mperf(const double *a, double *b);
  ARPREC_API void c_mpbessel(const double *a, double *b);
  ARPREC_API void c_mpbesselexp(const double *a, double *b);
  ARPREC_API void c_mpgamma(const double *a, double *b);

  /* assignment */
  ARPREC_API void c_mpeq(const double *a, double *b);
  ARPREC_API void c_mpeq_int(int a, double *b);
  ARPREC_API void c_mpeq_d(double a, double* b);
  ARPREC_API void c_mpeq_ji(int a, double *b);
  ARPREC_API void c_mpeq_zq(const double *a, double *b);
  ARPREC_API void c_mpeq_zx(double *r, double *i, double *b);
  ARPREC_API void c_mpeq_zz(const double *a, double *b);
  ARPREC_API void c_mpeq_str(const char *str, double *a);

  /* power */
  ARPREC_API void c_mppwr(const double *a, const double *b, double *c);
  ARPREC_API void c_mppwr_d(const double *a, double b, double *c);
  ARPREC_API void c_mppwr_qi(const double *a, int b, double *c);
  ARPREC_API void c_mppwr_jj(const double *a, const double *b, double *c);
  ARPREC_API void c_mppwr_ji(const double *a, int b, double *c);
  ARPREC_API void c_mppwr_zi(const double *a, int b, double *c);
  ARPREC_API void c_mppwr_zq(const double *a, const double *b, double *c);
  ARPREC_API void c_mppwr_zz(const double *a, const double *b, double *c);

  /* equality */
  ARPREC_API void c_mpcpr(const double *a, const double *b, int *c);
  ARPREC_API void c_mpcpr_i(const double *a, int b, int *c);
  ARPREC_API void c_mpcpr_d(const double *a, double b, int *c);
  ARPREC_API void c_mpcpr_z(const double *a, const double *b, int *c);

  /* less-than-or-equal-to */
  ARPREC_API void c_mplet(const double *a, const double *b, int *c);
  ARPREC_API void c_mplet_i(const double *a, int b, int *c);
  ARPREC_API void c_mplet_d(const double *a, double b, int *c);
  /* greater-than-or-equal-to */
  ARPREC_API void c_mpget(const double *a, const double *b, int *c);
  ARPREC_API void c_mpget_d(const double *a, double b, int *c);
  ARPREC_API void c_mpget_i(const double *a, int b, int *c);
  /* less-than */
  ARPREC_API void c_mpltt(const double *a, const double *b, int *c);
  ARPREC_API void c_mpltt_i(const double *a, int b, int *c);
  ARPREC_API void c_mpltt_d(const double *a, double b, int *c);
  /* greater-than */
  ARPREC_API void c_mpgtt(const double *a, const double *b, int *c);
  ARPREC_API void c_mpgtt_i(const double *a, int b, int *c);
  ARPREC_API void c_mpgtt_d(const double *a, double b, int *c);

  ARPREC_API void c_mpabs(const double *a, double *b);
  ARPREC_API void c_mpabs_z(const double *a, double *b);
  ARPREC_API void c_mparg(const double *a, double *b);

  ARPREC_API void c_mpacos(const double *a, double *b);
  ARPREC_API void c_mpasin(const double *a, double *b);
  ARPREC_API void c_mpatan(const double *a, double *b);
  ARPREC_API void c_mpatan2(const double *a, const double *b, double *c);
  ARPREC_API void c_mpcos(const double *a, double *b);
  ARPREC_API void c_mpcos_z(const double *a, double *b);
  ARPREC_API void c_mpcosh(const double *a, double *b);
  ARPREC_API void c_mpdble(const double *a, double *b);
  ARPREC_API void c_mpexp(const double *a, double *b);
  ARPREC_API void c_mpexp_z(const double *a, double *b);
  ARPREC_API void c_mpanint(const double *a, double *b);
  ARPREC_API void c_mpnint(const double *a, double *b);
  ARPREC_API void c_mplog(const double *a, double *b);
  ARPREC_API void c_mplog_z(const double *a, double *b);
  ARPREC_API void c_mplog10(const double *a, double *b);
  ARPREC_API void c_mpsin(const double *a, double *b);
  ARPREC_API void c_mpsin_z(const double *a, double *b);
  ARPREC_API void c_mpsinh(const double *a, double *b);
  ARPREC_API void c_mpnrt(const double *a, int *b, double *c);
  ARPREC_API void c_mpsqrt(const double *a, double *b);
  ARPREC_API void c_mpsqrt_z(const double *a, double *b);
  ARPREC_API void c_mptan(const double *a, double *b);
  ARPREC_API void c_mptanh(const double *a, double *b);
  ARPREC_API void c_mpmod(const double *a, const double *b, double *c);
  ARPREC_API void c_mpcssnf(const double *a, double *b, double *c);
  ARPREC_API void c_mpcsshf(const double *a, double *b, double *c);
  ARPREC_API void c_mprand(double *a);
  ARPREC_API void c_mpaint(const double *a, double *b);

  ARPREC_API void c_mpinp(const double *q);
  ARPREC_API void c_mpout(const double *q);
  ARPREC_API void c_mpwrite(const double *q, char *str, int *len);
  ARPREC_API void c_mpout_z(const double *q);
  ARPREC_API void c_mp_to_str(const double *a, char *str, int n_digits);

  ARPREC_API int c_mpinit(int nr_digits);
  ARPREC_API int c_mpinit_x(int nr_digits, const char *filename);
  ARPREC_API void c_mpgetpar(int pnum, int *val, int index);
  ARPREC_API void c_mpsetpar(int pnum, int val, int index);
  ARPREC_API void c_mpdotd(int *n, int *isa, double *a, int *isb, 
      const double *db, double *c);

  ARPREC_API void c_ovcheck(const double *a);
  ARPREC_API void c_mpinfr(const double *a, double *b, double *c);

  ARPREC_API void c_mpsetoutputprec(int num_digits);
  ARPREC_API int c_mpgetoutputprec(void);
  ARPREC_API void c_mpsetprec(int num_digits);
  ARPREC_API int c_mpgetprec(void);
  ARPREC_API void c_mpsetprecwords(int num_words);
  ARPREC_API int c_mpgetprecwords(void);

  ARPREC_API void c_mppi(double *a);

  ARPREC_API double *c_mpnew(void);
  ARPREC_API double *c_mpnew_x(int nr_words);
  ARPREC_API void c_mpfree(double *a);

#ifdef __cplusplus
}
#endif

#endif 

