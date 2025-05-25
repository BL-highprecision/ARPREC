/*
 * src/f_mp.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 */
#include "config.h"
#include <arprec/mp_real.h>
#include <arprec/c_mp.h>
#include <arprec/fpu.h>

#ifndef FC_FUNC_
#define FC_FUNC_(name, NAME)   name ## _
#endif

#define f_mpadd        FC_FUNC_(f_mpadd, F_MPADD)
#define f_mpadd_d      FC_FUNC_(f_mpadd_d, F_MPADD_D)
#define f_mpadd_ji     FC_FUNC_(f_mpadd_ji, F_MPADD_JI)
#define f_mpadd_jd     FC_FUNC_(f_mpadd_jd, F_MPADD_JD)
#define f_mpadd_zq     FC_FUNC_(f_mpadd_zq, F_MPADD_ZQ)
#define f_mpadd_zx     FC_FUNC_(f_mpadd_zx, F_MPADD_ZX)
#define f_mpadd_zz     FC_FUNC_(f_mpadd_zz, F_MPADD_ZZ)
#define f_mpsub        FC_FUNC_(f_mpsub, F_MPSUB)
#define f_mpsub_d      FC_FUNC_(f_mpsub_d, F_MPSUB_D)
#define f_mpsub_dq     FC_FUNC_(f_mpsub_dq, F_MPSUB_DQ)
#define f_mpsub_ji     FC_FUNC_(f_mpsub_ji, F_MPSUB_JI)
#define f_mpsub_ij     FC_FUNC_(f_mpsub_ij, F_MPSUB_IJ)
#define f_mpsub_jd     FC_FUNC_(f_mpsub_jd, F_MPSUB_JD)
#define f_mpsub_dj     FC_FUNC_(f_mpsub_dj, F_MPSUB_DJ)
#define f_mpsub_zq     FC_FUNC_(f_mpsub_zq, F_MPSUB_ZQ)
#define f_mpsub_qz     FC_FUNC_(f_mpsub_qz, F_MPSUB_QZ)
#define f_mpsub_zx     FC_FUNC_(f_mpsub_zx, F_MPSUB_ZX)
#define f_mpsub_xz     FC_FUNC_(f_mpsub_xz, F_MPSUB_XZ)
#define f_mpsub_zz     FC_FUNC_(f_mpsub_zz, F_MPSUB_ZZ)
#define f_mpneg_q      FC_FUNC_(f_mpneg_q, F_MPNEG_Q)
#define f_mpneg_z      FC_FUNC_(f_mpneg_z, F_MPNEG_Z)
#define f_mpmul        FC_FUNC_(f_mpmul, F_MPMUL)
#define f_mpmul_ji     FC_FUNC_(f_mpmul_ji, F_MPMUL_JI)
#define f_mpmul_qd     FC_FUNC_(f_mpmul_qd, F_MPMUL_QD)
#define f_mpmul_qi     FC_FUNC_(f_mpmul_qi, F_MPMUL_QI)
#define f_mpmul_zq     FC_FUNC_(f_mpmul_zq, F_MPMUL_ZQ)
#define f_mpmul_zz     FC_FUNC_(f_mpmul_zz, F_MPMUL_ZZ)
#define f_mpmul_zd     FC_FUNC_(f_mpmul_zd, F_MPMUL_ZD)
#define f_mpdiv        FC_FUNC_(f_mpdiv, F_MPDIV)
#define f_mpdiv_jj     FC_FUNC_(f_mpdiv_jj, F_MPDIV_JJ)
#define f_mpdiv_ji     FC_FUNC_(f_mpdiv_ji, F_MPDIV_JI)
#define f_mpdiv_ij     FC_FUNC_(f_mpdiv_ij, F_MPDIV_IJ)
#define f_mpdiv_qi     FC_FUNC_(f_mpdiv_qi, F_MPDIV_QI)
#define f_mpdiv_iq     FC_FUNC_(f_mpdiv_iq, F_MPDIV_IQ)
#define f_mpdiv_qd     FC_FUNC_(f_mpdiv_qd, F_MPDIV_QD)
#define f_mpdiv_dq     FC_FUNC_(f_mpdiv_dq, F_MPDIV_DQ)
#define f_mpdiv_zq     FC_FUNC_(f_mpdiv_zq, F_MPDIV_ZQ)
#define f_mpdiv_qz     FC_FUNC_(f_mpdiv_qz, F_MPDIV_QZ)
#define f_mpdiv_zz     FC_FUNC_(f_mpdiv_zz, F_MPDIV_ZZ)
#define f_mpdiv_zd     FC_FUNC_(f_mpdiv_zd, F_MPDIV_ZD)
#define f_mpdiv_dz     FC_FUNC_(f_mpdiv_dz, F_MPDIV_DZ)
#define f_mpdmc        FC_FUNC_(f_mpdmc, F_MPDMC)
#define f_mpmdc        FC_FUNC_(f_mpmdc, F_MPMDC)
#define f_mpeq         FC_FUNC_(f_mpeq, F_MPEQ)
#define f_mpeq_int     FC_FUNC_(f_mpeq_int, F_MPEQ_INT)
#define f_mpeq_d       FC_FUNC_(f_mpeq_d, F_MPEQ_D)
#define f_mpeq_ji      FC_FUNC_(f_mpeq_ji, F_MPEQ_JI)
#define f_mpeq_zq      FC_FUNC_(f_mpeq_zq, F_MPEQ_ZQ)
#define f_mpeq_zx      FC_FUNC_(f_mpeq_zx, F_MPEQ_ZX)
#define f_mpeq_zz      FC_FUNC_(f_mpeq_zz, F_MPEQ_ZZ)
#define f_mppwr        FC_FUNC_(f_mppwr, F_MPPWR)
#define f_mppwr_d      FC_FUNC_(f_mppwr_d, F_MPPWR_D)
#define f_mppwr_qi     FC_FUNC_(f_mppwr_qi, F_MPPWR_QI)
#define f_mppwr_jj     FC_FUNC_(f_mppwr_jj, F_MPPWR_JJ)
#define f_mppwr_ji     FC_FUNC_(f_mppwr_ji, F_MPPWR_JI)
#define f_mppwr_zi     FC_FUNC_(f_mppwr_zi, F_MPPWR_ZI)
#define f_mppwr_zq     FC_FUNC_(f_mppwr_zq, F_MPPWR_ZQ)
#define f_mppwr_zz     FC_FUNC_(f_mppwr_zz, F_MPPWR_ZZ)
#define f_mpcpr        FC_FUNC_(f_mpcpr, F_MPCPR)
#define f_mpcpr_i      FC_FUNC_(f_mpcpr_i, F_MPCPR_I)
#define f_mpcpr_d      FC_FUNC_(f_mpcpr_d, F_MPCPR_D)
#define f_mpcpr_z      FC_FUNC_(f_mpcpr_z, F_MPCPR_Z)
#define f_mplet        FC_FUNC_(f_mplet, F_MPLET)
#define f_mplet_i      FC_FUNC_(f_mplet_i, F_MPLET_I)
#define f_mplet_d      FC_FUNC_(f_mplet_d, F_MPLET_D)
#define f_mpget        FC_FUNC_(f_mpget, F_MPGET)
#define f_mpget_i      FC_FUNC_(f_mpget_i, F_MPGET_I)
#define f_mpget_d      FC_FUNC_(f_mpget_d, F_MPGET_D)
#define f_mpltt        FC_FUNC_(f_mpltt, F_MPLTT)
#define f_mpltt_i      FC_FUNC_(f_mpltt_i, F_MPLTT_I)
#define f_mpltt_d      FC_FUNC_(f_mpltt_d, F_MPLTT_D)
#define f_mpgtt        FC_FUNC_(f_mpgtt, F_MPGTT)
#define f_mpgtt_i      FC_FUNC_(f_mpgtt_i, F_MPGTT_I)
#define f_mpgtt_d      FC_FUNC_(f_mpgtt_d, F_MPGTT_D)

#define f_mpabs        FC_FUNC_(f_mpabs, F_MPABS)
#define f_mpabs_z      FC_FUNC_(f_mpabs_z, F_MPABS_Z)
#define f_mparg        FC_FUNC_(f_mparg, F_MPARG)
#define f_mpacos       FC_FUNC_(f_mpacos, F_MPACOS)
#define f_mpasin       FC_FUNC_(f_mpasin, F_MPASIN)
#define f_mpatan       FC_FUNC_(f_mpatan, F_MPATAN)
#define f_mpatan2      FC_FUNC_(f_mpatan2, F_MPATAN2)
#define f_mpcos        FC_FUNC_(f_mpcos, F_MPCOS)
#define f_mpcos_z      FC_FUNC_(f_mpcos_z, F_MPCOS_Z)
#define f_mpcosh       FC_FUNC_(f_mpcosh, F_MPCOSH)
#define f_mpdble       FC_FUNC_(f_mpdble, F_MPDBLE)
#define f_mpexp        FC_FUNC_(f_mpexp, F_MPEXP)
#define f_mpexp_z      FC_FUNC_(f_mpexp_z, F_MPEXP_Z)
#define f_mpaint       FC_FUNC_(f_mpaint, F_MPAINT)
#define f_mpnint       FC_FUNC_(f_mpnint, F_MPNINT)
#define f_mplog        FC_FUNC_(f_mplog, F_MPLOG)
#define f_mplog_z      FC_FUNC_(f_mplog_z, F_MPLOG_Z)
#define f_mplog10      FC_FUNC_(f_mplog10, F_MPLOG10)
#define f_mpsin        FC_FUNC_(f_mpsin, F_MPSIN)
#define f_mpsin_z      FC_FUNC_(f_mpsin_z, F_MPSIN_Z)
#define f_mpsinh       FC_FUNC_(f_mpsinh, F_MPSINH)
#define f_mpnrt        FC_FUNC_(f_mpnrt, F_MPNRT)
#define f_mpsqrt       FC_FUNC_(f_mpsqrt, F_MPSQRT)
#define f_mpsqrt_z     FC_FUNC_(f_mpsqrt_z, F_MPSQRT_Z)
#define f_mptan        FC_FUNC_(f_mptan, F_MPTAN)
#define f_mptanh       FC_FUNC_(f_mptanh, F_MPTANH)
#define f_mpmod        FC_FUNC_(f_mpmod, F_MPMOD)
#define f_mpcssnf      FC_FUNC_(f_mpcssnf, F_MPCSSNF)
#define f_mpcsshf      FC_FUNC_(f_mpcsshf, F_MPCSSHF)
#define f_mprand       FC_FUNC_(f_mprand, F_MPRAND)

#define f_mperfc       FC_FUNC_(f_mperfc, F_MPERFC)
#define f_mperf        FC_FUNC_(f_mperf, F_MPERF)
#define f_mpbessel     FC_FUNC_(f_mpbessel, F_MPBESSEL)
#define f_mpbesselexp  FC_FUNC_(f_mpbesselexp, F_MPBESSELEXP)
#define f_mpgamma      FC_FUNC_(f_mpgamma, F_MPGAMMA)

#define f_mpinp        FC_FUNC_(f_mpinp, F_MPINP)
#define f_mpout        FC_FUNC_(f_mpout, F_MPOUT)
#define f_mpout_z      FC_FUNC_(f_mpout_z, F_MPOUT_Z)
#define f_mpinit       FC_FUNC_(f_mpinit, F_MPINIT)
#define f_mpgetpar     FC_FUNC_(f_mpgetpar, F_MPGETPAR)
#define f_mpsetpar     FC_FUNC_(f_mpsetpar, F_MPSETPAR)
#define f_mpdotd       FC_FUNC_(f_mpdotd, F_MPDOTD)
#define f_ovcheck      FC_FUNC_(f_ovcheck, F_OVCHECK)
#define f_mpinfr       FC_FUNC_(f_mpinfr, F_MPINFR)
#define f_mp_to_str    FC_FUNC_(f_mp_to_str, F_MP_TO_STR)

#define f_mpsetoutputprec FC_FUNC_(f_mpsetoutputprec, F_MPSETOUTPUTPREC)
#define f_mpgetoutputprec FC_FUNC_(f_mpgetoutputprec, F_MPGETOUTPUTPREC)
#define f_mpsetprec       FC_FUNC_(f_mpsetprec, F_MPSETPREC)
#define f_mpgetprec       FC_FUNC_(f_mpgetprec, F_MPGETPREC)
#define f_mpsetprecwords  FC_FUNC_(f_mpsetprecwords, F_MPSETPRECWORDS)
#define f_mpgetprecwords  FC_FUNC_(f_mpgetprecwords, F_MPGETPRECWORDS)

#define f_fpu_fix_start  FC_FUNC_(f_fpu_fix_start, F_FPU_FIX_START)
#define f_fpu_fix_end   FC_FUNC_(f_fpu_fix_end,  F_FPU_FIX_STOP)

#ifdef __cplusplus
extern "C" {
#endif

static void f_mpinit_x(int n_digits, int compute_consts, int n_mantissa, int *n_words, 
    double *mpeps, double *_log2, double *_log10, double *_pi) {

  mp::mp_init(n_digits, NULL, static_cast<bool>(compute_consts));

  _log2[0] = _log10[0] = _pi[0] = mp::n_words;
  mp_real mpl022(_log2), mpl102(_log10), mppic2(_pi);
  mpl022 = mp_real::_log2;
  mpl102 = mp_real::_log10;
  mppic2 = mp_real::_pi;
  mpl022.toTempAndDestroy();
  mpl102.toTempAndDestroy();
  mppic2.toTempAndDestroy();

  mp::fmpwds5 = n_mantissa + 5;
  *n_words = mp::n_words;
  
  mpeps[0] = *n_words;
  mp_real mpeps2(mpeps);
  mpeps2 = mp_real::_eps;
  mpeps2.toTempAndDestroy();
}

void f_mpinit(const int *n_digits, const int *compute_consts, const int *n_mantissa, 
    int *n_words, double *mpeps, double *_log2, double *_log10, double *_pi) {
  f_mpinit_x(*n_digits, *compute_consts, *n_mantissa, n_words, mpeps, _log2, _log10, _pi);
}

void f_mpgetpar(const int *pnum, int *val, const int *index)
{
  c_mpgetpar(*pnum, val, *index);
}

void f_mpsetpar(const int *pnum, const int *val, const int *index)
{
  c_mpsetpar(*pnum, *val, *index);
}


/* add */
void f_mpadd(const double *a, const double *b, double *c) {
  c_mpadd(a, b, c);
}
void f_mpadd_d(const double *a, const double *b, double *c) {
  c_mpadd_d(a, *b, c);
}
void f_mpadd_ji(const double *a, const int *b, double *c) {
  c_mpadd_ji(a, *b, c);
}
void f_mpadd_jd(const double *a, const double *b, double *c) {
  c_mpadd_jd(a, *b, c);
}
void f_mpadd_zq(const double *a, const double *b, double *c) {
  c_mpadd_zq(a, b, c);
}
void f_mpadd_zx(const double *a, const double *br, 
                  const double *bi, double *c) {
  c_mpadd_zx(a, *br, *bi, c);
}
void f_mpadd_zz(const double *a, const double *b, double *c) {
  c_mpadd_zz(a, b, c);
}

/* sub */
void f_mpsub(const double *a, const double *b, double *c) {
  c_mpsub(a, b, c);
}

void f_mpsub_d(const double *a, const double *b, double *c) {
  c_mpsub_d(a, *b, c);
}
void f_mpsub_dq(const double *a, const double *b, double *c) {
  c_mpsub_dq(*a, b, c);
}
void f_mpsub_ji(const double *a, const int *b, double *c) {
  c_mpsub_ji(a, *b, c);
}
void f_mpsub_ij(const int *a, const double *b, double *c) {
  c_mpsub_ij(*a, b, c);
}
void f_mpsub_jd(const double *a, const double *b, double *c) {
  c_mpsub_jd(a, *b, c);
}
void f_mpsub_dj(const double *a, const double *b, double *c) {
  c_mpsub_dj(*a, b, c);
}
void f_mpsub_zq(const double *a, const double *b, double *c) {
  c_mpsub_zq(a, b, c);
}
void f_mpsub_qz(const double *a, const double *b, double *c) {
  c_mpsub_qz(a, b, c);
}
void f_mpsub_zx(const double *a, const double *br, 
                  const double *bi, double *c) {
  c_mpsub_zx(a, *br, *bi, c);
}
void f_mpsub_xz(const double *ar, const double *ai, 
                  const double *b, double *c) {
  c_mpsub_xz(*ar, *ai, b, c);
}
void f_mpsub_zz(const double *a, const double *b, double *c) {
  c_mpsub_zz(a, b, c);
}

/* negation */
void f_mpneg_q(const double *a, double *c) {
  c_mpneg_q(a, c);
}
void f_mpneg_z(const double *a, double *c) {
  c_mpneg_z(a, c);
}

/* mul */
void f_mpmul(const double *a, const double *b, double *c) {
  c_mpmul(a, b, c);
}
void f_mpmul_ji(const double *a, const int *b, double *c) {
  c_mpmul_ji(a, *b, c);
}
void f_mpmul_qi(const double *a, const int *b, double *c) {
  c_mpmul_qi(a, *b, c);
}
void f_mpmul_qd(const double *a, const double *b, double *c) {
  c_mpmul_qd(a, *b, c);
}
void f_mpmul_zq(const double *a, const double *b, double *c) {
  c_mpmul_zq(a, b, c);
}
void f_mpmul_zz(const double *a, const double *b, double *c) {
  c_mpmul_zz(a, b, c);
}
void f_mpmul_zd(const double *a, const double *b, double *c) {
  c_mpmul_zd(a, *b, c);
}

/* div */
void f_mpdiv(const double *a, const double *b, double *c) {
  c_mpdiv(a, b, c);
}
void f_mpdiv_jj(const double *a, const double *b, double *c) {
  c_mpdiv_jj(a, b, c);
}
void f_mpdiv_ji(const double *a, const int *b, double *c) {
  c_mpdiv_ji(a, *b, c);
}
void f_mpdiv_ij(const int *a, const double *b, double *c) {
  c_mpdiv_ij(*a, b, c);
}
void f_mpdiv_qi(const double *a, const int *b, double *c) {
  c_mpdiv_qi(a, *b, c);
}
void f_mpdiv_iq(const int *a, const double *b, double *c) {
  c_mpdiv_iq(*a, b, c);
}
void f_mpdiv_qd(const double *a, const double *b, double *c) {
  c_mpdiv_qd(a, *b, c);
}
void f_mpdiv_dq(const double *a, const double *b, double *c) {
  c_mpdiv_dq(*a, b, c);
}
void f_mpdiv_zq(const double *a, const double *b, double *c) {
  c_mpdiv_zq(a, b, c);
}
void f_mpdiv_qz(const double *a, const double *b, double *c) {
  c_mpdiv_qz(a, b, c);
}
void f_mpdiv_zz(const double *a, const double *b, double *c) {
  c_mpdiv_zz(a, b, c);
}
void f_mpdiv_zd(const double *a, const double *b, double *c) {
  c_mpdiv_zd(a, *b, c);
}
void f_mpdiv_dz(const double *a, const double *b, double *c) {
  c_mpdiv_dz(*a, b, c);
}

void f_mpdmc(const double *a, double *b) {
  c_mpdmc(*a, b);
}
void f_mpmdc(const double *a, double *b, int *n) {
  c_mpmdc(a, b, n);
}

/*special function*/
void f_mperfc(const double *a, double *b){
    c_mperfc(a,b);
}
void f_mperf(const double *a, double *b){
    c_mperf(a,b);
}
void f_mpbessel(const double *a, double *b){
    c_mpbessel(a,b);
}
void f_mpbesselexp(const double *a, double *b){
    c_mpbesselexp(a,b);
}
void f_mpgamma(const double *a, double *b){
    c_mpgamma(a,b);
}

/* assignment */
void f_mpeq(const double *a, double *b) {
  c_mpeq(a, b);
}
void f_mpeq_int(const int *a, double *b) {
  c_mpeq_int(*a, b);
}
void f_mpeq_d(const double *a, double *b) {
  c_mpeq_d(*a, b);
}
void f_mpeq_ji(const int *a, double *b) {
  c_mpeq_ji(*a, b);
}
void f_mpeq_zq(const double *a, double *b) {
  c_mpeq_zq(a, b);
}
void f_mpeq_zx(double *r, double *i, double *b) {
  c_mpeq_zx(r, i, b);
}
void f_mpeq_zz(const double *a, double *b) {
  c_mpeq_zz(a, b);
}


/* power */
void f_mppwr(const double *a, const double *b, double *c) {
  c_mppwr(a, b, c);
}
void f_mppwr_d(const double *a, const double *b, double *c) {
  c_mppwr_d(a, *b, c);
}
void f_mppwr_qi(const double *a, const int *b, double *c) {
  c_mppwr_qi(a, *b, c);
}
void f_mppwr_jj(const double *a, const double *b, double *c) {
  c_mppwr_jj(a, b, c);
}
void f_mppwr_ji(const double *a, const int *b, double *c) {
  c_mppwr_ji(a, *b, c);
}
void f_mppwr_zi(const double *a, const int *b, double *c) {
  c_mppwr_zi(a, *b, c);
}
void f_mppwr_zq(const double *a, const double *b, double *c) {
  c_mppwr_zq(a, b, c);
}

/* equality */
void f_mpcpr(const double *a, const double *b, int *c) {
  c_mpcpr(a, b,c);
}
void f_mpcpr_i(const double *a, const int *b, int *c) {
  c_mpcpr_i(a, *b, c);
}
void f_mpcpr_d(const double *a, const double *b, int *c) {
  c_mpcpr_d(a, *b, c);
}
void f_mpcpr_z(const double *a, const double *b, int *c) {
  c_mpcpr_z(a, b, c);
}

/* less-than-or-equal-to */
void f_mplet(const double *a, const double *b, int *c) {
  c_mplet(a, b, c);
}
void f_mplet_i(const double *a, const int *b, int *c) {
  c_mplet_i(a, *b, c);
}
void f_mplet_d(const double *a, const double *b, int *c) {
  c_mplet_d(a, *b, c);
}

/* greater-than-or-equal-to */
void f_mpget(const double *a, const double *b, int *c) {
  c_mpget(a, b, c);
}
void f_mpget_i(const double *a, const int *b, int *c) {
  c_mpget_i(a, *b, c);
}
void f_mpget_d(const double *a, const double *b, int *c) {
  c_mpget_d(a, *b, c);
}

/* less-than */
void f_mpltt(const double *a, const double *b, int *c) {
  c_mpltt(a, b, c);
}
void f_mpltt_i(const double *a, const int *b, int *c) {
  c_mpltt_i(a, *b, c);
}
void f_mpltt_d(const double *a, const double *b, int *c) {
  c_mpltt_d(a, *b, c);
}

/* greater-than */
void f_mpgtt(const double *a, const double *b, int *c) {
  c_mpgtt(a, b, c);
}
void f_mpgtt_i(const double *a, const int *b, int *c) {
  c_mpgtt_i(a, *b, c);
}
void f_mpgtt_d(const double *a, const double *b, int *c) {
  c_mpgtt_d(a, *b, c);
}

void f_mpabs(const double *a, double *b) {
  c_mpabs(a, b);
}
void f_mpabs_z(const double *a, double *b) {
  c_mpabs_z(a, b);
}

void f_mparg(const double *a, double *b) {
  c_mparg(a, b);
}

/* trigonometric functions */
void f_mpacos(const double *a, double *b) {
  c_mpacos(a, b);
}
void f_mpasin(const double *a, double *b) {
  c_mpasin(a, b);
}
void f_mpatan(const double *a, double *b) {
  c_mpatan(a, b);
}
void f_mpatan2(const double *a, const double *b, double *c) {
  c_mpatan2(a, b, c);
}
void f_mpcos(const double *a, double *b) {
  c_mpcos(a, b);
}
void f_mpcos_z(const double *a, double *b) {
  c_mpcos_z(a, b);
}
void f_mpdble(const double *a, double *b) {
  c_mpdble(a, b);
}
void f_mpcosh(const double *a, double *b) {
  c_mpcosh(a, b);
}
void f_mpexp(const double *a, double *b) {
  c_mpexp(a, b);
}
void f_mpexp_z(const double *a, double *b) {
  c_mpexp_z(a, b);
}

void f_mpaint(const double *a, double *b) {
  c_mpaint(a, b);
}
void f_mpnint(const double *a, double *b) {
  c_mpnint(a, b);
}

void f_mplog(const double *a, double *b) {
  c_mplog(a, b);
}
void f_mplog_z(const double *a, double *b) {
  c_mplog_z(a, b);
}
void f_mplog10(const double *a, double *b) {
  c_mplog10(a, b);
}
void f_mpsin(const double *a, double *b) {
  c_mpsin(a, b);
}
void f_mpsin_z(const double *a, double *b) {
  c_mpsin_z(a, b);
}
void f_mpsinh(const double *a, double *b) {
  c_mpsinh(a, b);
}
void f_mpnrt(const double *a, int *b, double *c) {
  c_mpnrt(a, b, c);
}
void f_mpsqrt(const double *a, double *b) {
  c_mpsqrt(a, b);
}
void f_mpsqrt_z(const double *a, double *b) {
  c_mpsqrt_z(a, b);
}

void f_mptan(const double *a, double *b) {
  c_mptan(a, b);
}
void f_mptanh(const double *a, double *b) {
  c_mptanh(a, b);
}

void f_mpmod(const double *a, const double *b, double *c) {
  c_mpmod(a, b, c);
}
void f_mpcsshf(const double *a, double *b, double *c) {
  c_mpcsshf(a, b, c);
}
void f_mpcssnf(const double *a, double *b, double *c) {
  c_mpcssnf(a, b, c);
}

void f_mprand(double *a) {
  c_mprand(a);
}

void f_mp_to_str(const double *a, char *s, int n_digits) {
  c_mp_to_str(a, s, n_digits);
}

void f_ovcheck(const double *a) {
  c_ovcheck(a);
}

void f_mpinfr(const double *a, double *b, double *c) {
  c_mpinfr(a, b, c);
}

/* Input */
void f_mpinp(const double *q) {
  c_mpinp(q);
}

/* Output */
void f_mpout(const double *q, char *c, int *l) {
  c_mpwrite(q, c, l);
}

void f_mpout_z(const double *q) {
  c_mpout_z(q);
}


void f_mpdotd(int *n, int *isa, double *a, int *isb, const double *db,
	      double *c) {
  c_mpdotd(n, isa, a, isb, db, c);
}

void f_fpu_fix_start(unsigned int *old_cw) {
  fpu_fix_start(old_cw);
}

void f_fpu_fix_end(unsigned int *old_cw) {
  fpu_fix_end(old_cw);
}

void f_mpsetoutputprec(int *num_digits) {
  c_mpsetoutputprec(*num_digits);
}

void f_mpgetoutputprec(int *num_digits) {
  *num_digits = c_mpgetoutputprec();
}

void f_mpsetprec(int *num_digits) {
  c_mpsetprec(*num_digits);
}

void f_mpgetprec(int *num_digits) {
  *num_digits = c_mpgetprec();
}

void f_mpsetprecwords(int *num_words) {
  c_mpsetprecwords(*num_words);
}

void f_mpgetprecwords(int *num_words) {
  *num_words = c_mpgetprecwords();
}

#ifdef __cplusplus
}
#endif

