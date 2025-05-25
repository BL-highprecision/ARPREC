#include <arprec/mp_real.h>

mp_real::mp_real(const char *s, unsigned int sz) {
  allocate(sz);
  if (s) {
    if (!read(s))
      std::cerr << "mp_real initialization from C string failed." << std::endl;
  }
}

mp_real::mp_real(const std::string &s, unsigned int sz) {
  allocate(sz);
  if (sz) {
    if (!read(s))
      std::cerr << "mp_real initialization from std::string failed." << std::endl;
  }
}

mp_real_temp pow(const mp_real& a, const mp_real& b) {
  mp_real t1, t2;
  int prec_words = mp::prec_words;
  mp_real::mplogx(a, mp_real::_pi, mp_real::_log2, t1, prec_words);
  mp_real::mpmulx(t1, b, t2, prec_words);
  mp_real::mpexpx(t2, mp_real::_pi, mp_real::_log2, t1);
  return t1.toTempAndDestroy();
}

mp_real_temp pow(const mp_real& a, double b) {
  mp_real t1, t2;
  int prec_words = mp::prec_words;
  mp_real::mplogx(a, mp_real::_pi, mp_real::_log2, t1, prec_words);
  mp_real::mpmuld(t1, b, 0, t2, prec_words);
  mp_real::mpexpx(t2, mp_real::_pi, mp_real::_log2, t1);
  return t1.toTempAndDestroy();
}

mp_real_temp pow(const mp_real& a, int n) {
  mp_real c;
  mp_real::mpnpwx(a, n, c, mp::prec_words);
  return c.toTempAndDestroy();
}

mp_real_temp acos(const mp_real& a) {
  mp_real c, temp, temp2,  f(1.0, 6);
  int prec_words = mp::prec_words;
  // compute b = sqrt(1- a^2);
  mp_real::mpmulx(a, a, temp, prec_words);
  mp_real::mpsub(f, temp, temp2, prec_words);
  mp_real::mpsqrtx(temp2, temp, prec_words);
  mp_real::mpangx(a, temp, mp_real::_pi, c);
  return c.toTempAndDestroy();
}

mp_real_temp aint(const mp_real& a) {
  mp_real ret, junk(0.0, 0);
  mp_real::mpinfr(a, ret, junk, mp::prec_words, 0);
  return ret.toTempAndDestroy();
}

mp_real_temp anint(const mp_real& a) {
  mp_real ret;
  mp_real::mpnint(a, ret, mp::prec_words);
  return ret.toTempAndDestroy();
}

mp_real_temp asin(const mp_real& a) {
  mp_real c, temp, temp2,  f(1.0, 6);
  int prec_words = mp::prec_words;
  // compute b = sqrt(1- a^2);
  mp_real::mpmulx(a, a, temp, prec_words);
  mp_real::mpsub(f, temp, temp2, prec_words);
  mp_real::mpsqrtx(temp2, temp, prec_words);
  mp_real::mpangx(temp, a, mp_real::_pi, c);
  return c.toTempAndDestroy();
}

mp_real_temp atan(const mp_real& a) {
  mp_real c;
  mp_real::mpangx(mp_real(1.0), a, mp_real::_pi, c);
  return c.toTempAndDestroy();
}

mp_real_temp atan2(const mp_real& y, const mp_real& x) {
  mp_real c;
  //XSL mp_real::mpangx(y, x, mp_real::_pi, c);
  mp_real::mpangx(x, y, mp_real::_pi, c);
  return c.toTempAndDestroy();
}

mp_real_temp cos(const mp_real& a) {
  mp_real c, junk;
  mp_real::mpcssx(a, mp_real::_pi, c, junk);
  return c.toTempAndDestroy();
}

mp_real_temp cosh(const mp_real& a) {
  mp_real c, junk;
  mp_real::mpcshx(a, mp_real::_pi, mp_real::_log2, c, junk);
  return c.toTempAndDestroy();
}

mp_real_temp exp(const mp_real& a) {
  mp_real c;
  mp_real::mpexpx(a, mp_real::_pi, mp_real::_log2, c);
  return c.toTempAndDestroy();
}

mp_real_temp log(const mp_real& a) {
  mp_real c;
  mp_real::mplogx(a, mp_real::_pi, mp_real::_log2, c, mp::prec_words);
  return c.toTempAndDestroy();
}

mp_real_temp log10(const mp_real& a) {
  mp_real c;
  int prec_words = mp::prec_words;
  mp_real::mplogx(a, mp_real::_pi, mp_real::_log2, c, prec_words);
  mp_real::mpdivx(c, mp_real::_log10, c, prec_words);
  return c.toTempAndDestroy();
}

mp_real_temp mp_rand() {
  mp_real a;
  mp_real::mprand(a);
  return a.toTempAndDestroy();
}

void mpcsshf(const mp_real& a, mp_real& b, mp_real& c) {
  mp_real::mpcshx(a, mp_real::_pi, mp_real::_log2, b, c);
  return;
}

void mpcssnf(const mp_real& a, mp_real& b, mp_real& c) {
  mp_real::mpcssx(a, mp_real::_pi, b, c);
  return;
}

mp_real_temp sin(const mp_real& a) {
  mp_real c, junk;
  mp_real::mpcssx(a, mp_real::_pi, junk, c);
  return c.toTempAndDestroy();
}

mp_real_temp sinh(const mp_real& a) {
  mp_real c, junk;
  mp_real::mpcshx(a, mp_real::_pi, mp_real::_log2, junk, c);
  return c.toTempAndDestroy();
}

mp_real_temp sqrt(const mp_real& a) {
  mp_real ret;
  mp_real::mpsqrtx(a, ret, mp::prec_words);
  return ret.toTempAndDestroy();
}

mp_real_temp sqr(const mp_real& a) {
  mp_real ret;
  mp_real::mpsqx(a, ret, mp::prec_words);
  return ret.toTempAndDestroy();
}

mp_real_temp tan(const mp_real& a) {
  mp_real c, d;
  mp_real::mpcssx(a, mp_real::_pi, d, c);
  return c / d;
}

mp_real_temp tanh(const mp_real& a) {
  mp_real c, d;
  mp_real::mpcshx(a, mp_real::_pi, mp_real::_log2, d, c);
  return c / d;
}

// Some compilers (IBM AIX xlC 6 in particular) needs this
// Other compilers automatically convert mp_real_tmp to mp_real.
std::ostream &operator<<(std::ostream &s, const mp_real_temp &x) {
  return (s << mp_real(x));
}

void mp_real::allocate(unsigned int sz) {
  if (sz == 0) { 
    mpr = NULL; 
    alloc = false;
  } else {
    alloc = true;
    mpr = new double[sz];
    mpr[0] = static_cast<double>(sz);
    mpr[1] = 0.0;
#if (ARPREC_DEBUG)
    for (unsigned int i = 2; i < sz; i++)
      mpr[i] = mp::_d_nan;
    VALGRIND_MAKE_MEM_UNDEFINED(mpr + 2, (sz-2)*sizeof(double));
#endif
  }
}

