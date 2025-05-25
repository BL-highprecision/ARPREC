#include <arprec/mp_real.h>
#include "small_inline.h"

mp_real_temp fmod(const mp_real& a, const mp_real &b)
{
  mp_real c, d, e;
  int nws = c.prec_words;
  int prec_words = mp::prec_words;

  c.prec_words = std::min(c.prec_words, int(a[2]) - int(b[2]) + 3);
  c.prec_words = std::max(c.prec_words, 1);
  mp_real::mpdivx(a, b, c, prec_words);
  c.prec_words = nws;
  mp_real::mpinfr(c, d, e, prec_words);
  mp_real::mpmulx(d, b, e, prec_words);
  mp_real::mpsub(a, e, c, prec_words);  /// THIS LINE WAS CHANGED

  //now make sure that division was correct
  //it is possible that it was slightly wrong.
  double a2 = a[1];
  double c2 = c[1];
  double b2 = b[1];
  if(a2 > 0.0) {
    if(c2 < 0.0) {
      mp_real::mpadd(c, (b>0.0) ? b : mp_real(-b), c, prec_words);
    } else if(b2 > 0.0) {
      if(c >= b) {
        mp_real::mpsub(c, b, c, prec_words);
      }
    } else if(b2 < 0.0) {
      c[1] = -c[1];
      if(c <= b) {
        mp_real::mpsub(c, b, c, prec_words);
      }
      c[1] = -c[1];
    }
  } else {
    if(c2 > 0.0) {
      mp_real::mpsub(c, (b>0.0) ? b : mp_real(-b), c, prec_words);
    } else if(b2 < 0.0) {
      if(c <= b) {
        mp_real::mpsub(c, b, c, prec_words);
      }
    } else if(b2 > 0.0) {
      c[1] = -c[1];
      if(c >= b) {
        mp_real::mpsub(c, b, c, prec_words);
      }
      c[1] = -c[1];
    }
  }
  return c.toTempAndDestroy();
}
