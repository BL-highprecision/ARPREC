/*
 * src/mpreal2.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 */
#include <arprec/mp_real.h>
#include "small_inline.h"
#include <cstdio>

using std::cerr;
using std::endl;

void mp_real::mpdotd(int n, int isa, const mp_real a[],
                     int isb, const double db[], mp_real& c)
{
  // This routine computes the dot product of the MP vector A with the DP
  // vector DB, returning the MP result in C.  This routine is used in the
  // author's customized PSLQ routine, resulting in substantial speedup.
  // The length of both the A and DB vectors is N, and ISA and ISB are the 
  // skip distances between successive elements of A and DB, measured in 
  // mp_real words, and in DP words, respectively. The DP values in DB must
  // be whole numbers, so for example they cannot be larger than 2^53 on
  // IEEE systems.  Debug output begins with debug_level = 8.
  //
  //mp_real s1((int)(prec_words+6));
  double *s1;
  double *d1 = new double[prec_words+6], *d2 = new double[prec_words+6];
  double dt0, dt1, dt2, t[2];
  int i, ixd, ish, k, ka, kb, kz, m1, m2, m3, m4, m5, nd, nsh, n1;
  const int mpnpr4 = std::max(0 , (mp_real::mpnpr >> 2) - 1);
  int prec_words = mp::prec_words;

  // Set DMAX to the largest DP whole number that can be represented exactly:
  // 2.d0 ** 53 on IEEE systems, or 2.d0 ** 48 on Cray non-IEEE systems.
#define MPDOT_CHECK_FOR_ERROR 
#ifdef MPDOTD_CHECK_FOR_ERROR
  double dmax = 9007199254740992.0; // 2^53
#endif
  //double two30 = 1073741824.0;      // 2^30

  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(c);
    return;
  }

  if (debug_level >= 8) {
    printf("mpdotd input: n = %d, isa = %d, isb = %d\n", n, isa, isb);
    for (k = 0; k < n; ++k) {
      ka = k * isa;
      kb = k * isb;
      mpmdc(a[ka], dt1, n1, prec_words);
      dt1 = ldexp(dt1, n1);
      printf("k = %d, dt1 = %25.16e, db = %25.16e\n", k, dt1, db[kb]);
    }
  }

  double two48 = 281474976710656.0; // 2^48
  for (i = 0; i < n; ++i) {
    kb = i * isb;
    if ( fabs(db[kb]) > two48 ) {  // use conventional loop
      // cout << "naive dot loop is called\n";
      mp_real s1, s2;
      zero(c);
      for (k = 0; k < n; ++k) {
        ka = k * isa;
        kb = k * isb;
        mpmuld (a[ka], db[kb], 0, s1, prec_words);
        mpadd (s1, c, s2, prec_words);
        mpeq (s2, c, prec_words);
      }
      return;
    }
  }
  
  for (i = 1; i <= 4; ++i) d2[i] = 0.;
  for (i = 0; i < prec_words+6; ++i) d1[i] = 0.;

  // ND is the length of D1, IXD is the exponent as in ordinary arprec format.
  // In the code below ND + 1 mantissa words are maintained whenever possible.
  nd = 0;
  ixd = 0;

  // Loop over the n input data pairs.

  for (k = 0; k < n; ++k) {
    ka = k * isa;
    kb = k * isb;

    s1 = a[ka].mpr;
    
    int na = std::min (std::abs(static_cast<int>(s1[1])), prec_words); // number of words in a[k]
    dt0 = db[kb];

    if (na == 0 || dt0 == 0.0) continue; // This pair is zero.

    // Check to make sure the input DP value satisfies the requirements.
#ifdef MPDOTD_CHECK_FOR_ERROR
    /* Check not used */
    if (std::abs (dt0) >= dmax || fmod (dt0, 1.0) != 0.) {
      if (MPKER[73] != 0) {
        printf((char*)"mpdotd: improper dp value: k = %d, dt0 = %25.15e\n", k, dt0);
        if (MPKER[73] == 2) mpabrt();
      }
      return;
      }
#endif


    // Save the two initial words of A.

    int ia1 = static_cast<int>(s1[1]);
    int ia2 = static_cast<int>(s1[2]);
    if (ia1 < 0) dt0 = - dt0;

    // Split the input DP value into high-order (<= 53-mpnbt bits) and 
    // low-order (<= mpnbt bits) values.
    dt1 = int (mprdx * dt0);
    dt2 = dt0 - mpbdx * dt1;

    if ( dt1 == 0.0 ) {
      // Only the low-order part of the DP value is nonzero.
      ish = ia2 - ixd;
      if (nd == 0) ish = 0;

      if (ish >= 0) {
        // The product a[k] * db[k] has greater exponent than the cumulative
        // sum. Thus the cumulative sum must be shifted to the right by
        // ish words.
        m1 = std::min (na, ish);
        m2 = std::min (na, nd + ish);
        m3 = na;
        m4 = std::min (std::max (na, ish), prec_words + 1);
        m5 = std::min (std::max (na, nd + ish), prec_words + 1);
        d2[1] = 0.;
        d2[2] = 0.;

        for (i = FST_M; i < m1 + FST_M; ++i) {
          t[0] = mp_two_prod(dt2, s1[i], t[1]);
          d2[i-1] += t[0];
          d2[i] = t[1];
        }

        for (i = m1 + FST_M; i < m2 + FST_M; ++i) {
          t[0] = mp_two_prod(dt2, s1[i], t[1]);
          d2[i-1] += t[0];
          d2[i] = t[1] + d1[i-ish];
        }

        for (i = m2 + FST_M; i < m3 + FST_M; ++i) {
          t[0] = mp_two_prod(dt2, s1[i], t[1]);
          d2[i-1] += t[0];
          d2[i] = t[1];
        }

        for (i = m3 + FST_M; i < m4 + FST_M; ++i)
          d2[i] = 0.0;
      
        for (i = m4 + FST_M; i < m5 + FST_M; ++i)
          d2[i] = d1[i-ish];

        // Copy d2 back to d1.
        if (d2[2] != 0) {
          ish = 1;
          ++m5;
          ixd = ia2 + 1;
        } else {
          ish = 0;
          ixd = ia2;
        }
        d1[2] = 0.0;
        for (i = 2; i < m5 + FST_M; ++i) d1[i+ish] = d2[i];

        nd = m5;
        d1[nd+3] = 0.0;
        d1[nd+4] = 0.0;

      } else {
        // The product a(k) * db(k) has smaller exponent than the cumulative
        // sum. Thus the product must be shifted to the right by -ish words.
        nsh = -ish;
        m1 = std::min (nd, nsh);
        m2 = std::min (nd, na + nsh);
        m3 = nd;
        m4 = std::min (std::max (nd, nsh), prec_words + 1);
        m5 = std::min (std::max (nd, na + nsh), prec_words + 1);
        
        for (i = m1 + FST_M; i < m2 + FST_M; ++i) {
          t[0] = mp_two_prod(dt2, s1[i-nsh], t[1]);
          d1[i-1] += t[0];
          d1[i] += t[1];
        }
        
        for (i = m3 + FST_M; i < m4 + FST_M; ++i)
          d1[i] = 0.0;

        for (i = m4 + FST_M; i < m5 + FST_M; ++i) {
          t[0] = mp_two_prod(dt2, s1[i-nsh], t[1]);
          d1[i-1] += t[0];
          d1[i] = t[1];
        }

        nd = m5;
        d1[nd+3] = 0.0;
        d1[nd+4] = 0.0;
      }

    } else {
      // All two parts of the input DP value are nonzero.
      ish = ia2 + 1 - ixd; // The product's exponent is incremented by 1.
      if (nd == 0) ish = 0;
      double s11 = s1[1];
      double s12 = s1[2];
      s1[1] = s1[2] = 0.0;

      if (ish >= 0) {
        // The product a[k] * db[k] has greater exponent that the cumulative
        // sum. Thus the cumulative sum must be shifted to the right by
        // ish words.
        m1 = std::min (na + 1, ish);
        m2 = std::min (na + 1, nd + ish);
        m3 = na + 1;
        m4 = std::min (std::max (na + 1, ish), prec_words + 1);
        m5 = std::min (std::max (na + 1, nd + ish), prec_words + 1);
        d2[1] = 0.0;
        d2[2] = 0.0;
        
        for (i = FST_M; i < m1 + FST_M; ++i) {
          t[0] = mp_two_prod(dt1, s1[i], t[1]);   // high-order part
          d2[i-1] += t[0];
          d2[i] = t[1];
          t[0] = mp_two_prod(dt2, s1[i-1], t[1]); // low-order part
          d2[i-1] += t[0];
          d2[i] += t[1];
        }

        for (i = m1 + FST_M; i < m2 + FST_M; ++i) {
          t[0] = mp_two_prod(dt1, s1[i], t[1]);   // high-order part
          d2[i-1] += t[0];
          d2[i] = t[1] + d1[i-ish];
          t[0] = mp_two_prod(dt2, s1[i-1], t[1]); // low-order part
          d2[i-1] += t[0];
          d2[i] += t[1];
        }

        for (i = m2 + FST_M; i < m3 + FST_M; ++i) {
          t[0] = mp_two_prod(dt1, s1[i], t[1]);   // high-order part
          d2[i-1] += t[0];
          d2[i] = t[1];
          t[0] = mp_two_prod(dt2, s1[i-1], t[1]); // low-order part
          d2[i-1] += t[0];
          d2[i] += t[1];
        }

        for (i = m3 + FST_M; i < m4 + FST_M; ++i)
          d2[i] = 0.0;
      
        for (i = m4 + FST_M; i < m5 + FST_M; ++i)
          d2[i] = d1[i-ish];

        // Copy d2 back to d1.
        if (d2[2] != 0) {
          ish = 1;
          ++m5;
          ixd = ia2 + 2;
        } else {
          ish = 0;
          ixd = ia2 + 1;
        }
        d1[2] = 0.0;
        for (i = 2; i < m5 + FST_M; ++i) d1[i+ish] = d2[i];

        nd = m5;
        d1[nd+3] = 0.0;
        d1[nd+4] = 0.0;

      } else {
        // The product a[k] * db[k] has smaller exponent that the cumulative
        // sum. Thus the product must be shifted to the right by -ish words.
        nsh = -ish;
        m1 = std::min (nd, nsh);
        m2 = std::min (nd, na + 1 + nsh);
        m3 = nd;
        m4 = std::min (std::max (nd, nsh), prec_words + 1);
        m5 = std::min (std::max (nd, na + 1 + nsh), prec_words + 1);

        for (i = m1 + FST_M; i < m2 + FST_M; ++i) {
          t[0] = mp_two_prod(dt1, s1[i-nsh], t[1]); // high-order part
          d1[i-1] += t[0];
          d1[i] += t[1];
          t[0] = mp_two_prod(dt2, s1[i-1-nsh], t[1]); // low-order part
          d1[i-1] += t[0];
          d1[i] += t[1];
        }

        for (i = m3 + FST_M; i < m4 + FST_M; ++i)
          d1[i] = 0.0;
        
        for (i = m4 + FST_M; i < m5 + FST_M; ++i) {
          t[0] = mp_two_prod(dt1, s1[i-nsh], t[1]); // high-order part
          d1[i-1] += t[0];
          d1[i] = t[1];
          t[0] = mp_two_prod(dt2, s1[i-1-nsh], t[1]); // low-order part
          d1[i-1] += t[0];
          d1[i] += t[1];
        }

        nd = m5;
        d1[nd+3] = 0.0;
        d1[nd+4] = 0.0;
      }
      s1[1] = s11;
      s1[2] = s12;
    }
    
    if (nd == 0) continue;

    // Release carry to avoid overflowing the exact integer capacity
    // (2^52-1) of a floating poitn word in D1.
    // Results may be negative, but that is not a problem -- these will be
    // fixed in the final call to mpnorm.
    if(!((k-1) & mpnpr4)) {
      dt2 = 0.0; // carry
      for (i = nd+2; i >= FST_M; --i) {
        dt1 = dt2 + d1[i];
        dt2 = int (dt1 * mprdx);   // carry < 2^(48 - mpnbt)
        d1[i] = dt1 - dt2 * mpbdx; // remainder of dt1 * 2^(-mpnbt)
      }
      d1[2] += dt2;
    }
    // If d1[2] is nonzero due to carry release, shift result to right.

    if (d1[2] != 0.0) {
      ish = 1;
      ++ixd;
      nd = std::min (nd + 1, prec_words + 1);
    } else {
      ish = 0;
    }

    if (ish != 0) {
      for (i = nd+2; i >= FST_M; --i) d1[i] = d1[i-ish];
      d1[1] = 0.0;
      d1[2] = 0.0;
    }
    d1[nd+3] = 0.0;
    d1[nd+4] = 0.0;
 
    // Check to see if there are leading zeros.
    for (i = FST_M; i < nd + FST_M; ++i)
      if (d1[i] != 0.0) break;

    if ( i == nd + FST_M) {
      // The cumulative sum is now zero.
      nd = 0;
      ixd = 0;
      d1[1] = 0.0;
      d1[2] = 0.0;
    } else {
      kz = i - FST_M;
      if (kz > 0) {
        // Leading zeroes -- shift cumulative sum to left.
        for (i = FST_M; i < nd - kz + FST_M; ++i)
          d1[i] = d1[i+kz];

        nd = nd - kz;
        ixd = ixd - kz;
        d1[nd+3] = 0.0;
        d1[nd+4] = 0.0;
      }
    }

  } // end for k = ...


  // Call mpnorm to fix up result and store in c.

  d1[1] = nd;
  d1[2] = ixd;

  mpnorm(d1, c, prec_words);

  delete [] d1;
  delete [] d2;
  
  if (debug_level >= 8) {
    mpmdc(c, dt1, n1, prec_words);
    dt1 = ldexp(dt1, n1);
    cerr << "mpdotd output: dt1 = " << dt1 << endl;
    print_mpreal("MPDOTD O: c ", c);
  }
}

