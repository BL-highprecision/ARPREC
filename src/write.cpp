/*
 * src/write.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2006
 *
 */
#include <string>
#include <sstream>
#include <iomanip>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::ostream;
using std::istream;
using std::cerr;
using std::endl;
using std::string;
using std::ios_base;

/* Custom streambuf class to limit each line to width characters. */
template <class Ch, class Tr = std::char_traits<Ch> >
class fixed_col_streambuf : public std::basic_streambuf<Ch, Tr> {
  std::basic_streambuf<Ch, Tr> *_sb;
  int _width;
  int _col;

  typedef std::basic_streambuf<Ch, Tr> base_type;
  typedef typename base_type::int_type int_type;
public:
  fixed_col_streambuf(base_type *sb, int width = 72) :
    _sb(sb), _width(width), _col(0) { }

  virtual int_type overflow(int_type ch = Tr::eof()) { 
    if (ch == Tr::eof()) return 0;

    if (_width > 0 && ch != '\n') {
      if (++_col > _width) {
        _col = 1;
        _sb->sputc('\n');
      }
    }
    return _sb->sputc(ch);
  }
};

template <class Ch, class Tr = std::char_traits<Ch> >
class fixed_col_ostream : public std::basic_ostream<Ch, Tr> {
  fixed_col_streambuf<Ch, Tr> *buf;
  typedef std::basic_ostream<Ch, Tr> base_type;
public:
  fixed_col_ostream(base_type &os, int width = 0) : base_type(new fixed_col_streambuf<Ch, Tr>(os.rdbuf(), width)) { }
  fixed_col_ostream(std::streambuf *buf, int width = 0) : base_type(new fixed_col_streambuf<Ch, Tr>(buf)) { }
  ~fixed_col_ostream() { delete this->rdbuf(); }
};

string mp_real::to_string(int precision, int width, int expn_width, 
    ios_base::fmtflags fmt, bool showpos, bool uppercase, char fill) const {
  std::ostringstream os;
  bool fixed = (fmt & ios_base::fixed) != 0;
  bool scientific = (fmt & ios_base::scientific) != 0;
  bool custom = (!fixed && !scientific);

  if (custom) {
    // Output in the custom   10 ^ expn x digits    format.
    std::ostringstream os;
    char *t = new char[precision + 1];
    int len, e;
    len = to_digits(t, e, precision);
    os << "10 ^ ";
    if (expn_width > 0) os << std::setw(expn_width);
    os << e << " x ";
    if (*this < 0.0) os << '-'; 
    else if (showpos) os << '+';
    os << t[0];
    if (len > 1)
      os << '.';
    os << (t+1);
    delete [] t;

    string s = os.str();

    /* Fill in the blanks */
    len = static_cast<int>(s.length());
    if (len < width) {
      int delta = width - len;
      if (fmt & ios_base::left) {
        s.append(delta, fill);
      } else {
        s.insert(static_cast<string::size_type>(0), delta, fill);
      }
    }

    return s;
  }

  bool sgn = true;
  int i, e = 0;

  if (*this < 0.0)
    os << '-';
  else if (showpos)
    os << '+';
  else
    sgn = false;

  if (*this == 0.0) {
    /* Zero case */
    os << '0';
    if (precision > 0) {
      os << '.';
      for (i = 0; i < precision; i++) os << '0';
    }
  } else {
    /* Non-zero case */
    int off = (fixed ? (1 + static_cast<int>((floor(dble(log10(abs(*this))))))) : 1);
    int d = precision + off;

    if (fixed && d <= 0) {
      os << '0';
      if (precision > 0) {
        os << '.';
        for (i = 0; i < precision; i++) os << '0';
      }
    } else {
      char *t = new char[d+1];
      int j;

      int len = to_digits(t, e, d);
      for (i = len; i <= d; i++) { t[i] = '0'; }

      if (fixed) {
        if (off > 0) {
          for (i = 0; i < off; i++) os << t[i];
          if (precision > 0) {
            os << '.';
            for (j = 0; j < precision; j++, i++) os << t[i];
          }
        } else {
          os << "0.";
          if (off < 0) for (i = 0; i < -off; i++) os << '0';
          for (i = 0; i < d; i++) os << t[i];
        }
      } else {
        os << t[0];
        if (precision > 0) os << '.';
        for (i = 1; i <= precision; i++) os << t[i];

        delete [] t;
      }
    }
  }

  if (!fixed) {
    /* Fill in exponent part */
    os << (uppercase ? 'E' : 'e');
    os << ((e >= 0) ? '+' : '-');
    os << std::setfill('0') << std::setw(2) << std::abs(e);
  }

  string s = os.str();

  /* Fill in the blanks */
  int len = static_cast<int>(s.length());
  if (len < width) {
    int delta = width - len;
    if (fmt & ios_base::internal) {
      if (sgn)
        s.insert(static_cast<string::size_type>(1), delta, fill);
      else
        s.insert(static_cast<string::size_type>(0), delta, fill);
    } else if (fmt & ios_base::left) {
      s.append(delta, fill);
    } else {
      s.insert(static_cast<string::size_type>(0), delta, fill);
    }
  }

  return s;
}

/* Writes a mp_real number to the output stream.
 * Number is written out in the form
 *   10 ^ expn x digits                             */
bool mp_real::write(ostream &s, int precision, int width, int expn_width, 
    ios_base::fmtflags fmt, int n_cols, 
    bool showpos, bool uppercase, char fill) const {
  if (n_cols > 0) {
    fixed_col_ostream<char> os(s, n_cols);
    return write(os, precision, width, expn_width, 
        fmt, 0, showpos, uppercase, fill);
  }

  string str = to_string(precision, width, expn_width, 
      fmt, showpos, uppercase, fill);
  return (s << str).good();
}

