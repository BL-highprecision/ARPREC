/*
 * src/read.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2006
 *
 */
#include <sstream>
#include <string>
#include <iomanip>
#include <cctype>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::string;
using std::istream;

static string read_digits(istream &s, bool allow_sign = true, bool skip_space = false) {
  string str;
  char ch;
  s >> std::ws;
  if (allow_sign) {
    s >> ch;
    if (ch == '-') {
      str = "-";
      if (skip_space) s >> std::ws;
    } else if (ch != '+')
      s.putback(ch);
  }

  for (;;) {
    if (skip_space) s >> std::ws;

    s >> ch;
    if (!s) break;

    if (isdigit(ch))
      str += ch;
    else {
      s.putback(ch);
      break;
    }
  }
  if (str == "-") str = "";
  return str;
}

/* Reads a mp_real from input stream.  Two possible format is as follows:
 *
 *   10 ^ expn x digits.digits
 *   digits.digits e expn
 *
 * In both cases exponent portion is optional.  The decimal point is also
 * optional, and cases does not matter.  For negative numbers, a minus sign
 * should appear in front of the first digit.
 *
 * White space in between digits is ignored, so multi-line numbers can be
 * entered.  To end the series of digits, one can use comma (,) or semicolon (;).
 *
 * example:  10 ^ 0 x 3.141592653589793238462643383279502884
 *                      197169399375105820974944592307816406
 *                      286208998628034825342117067982148086           */
bool mp_real::read(istream &s) {
  char ch;
  bool expn_before = true;
  string str;
  string expn_str;
  string digit_str;

  istream::fmtflags old_flags = s.flags();
  s >> std::noskipws;
  s >> std::ws;

  ch = s.peek();
  if (ch == '+' || ch == '-') expn_before = false;
  str = read_digits(s);

  if (expn_before && str != "10") expn_before = false;
  
  s >> std::ws;
  s >> ch;
  if (expn_before && ch == '^') {
    str = "";
    expn_str = read_digits(s);
    s >> std::ws;
    s >> ch;
    if (ch == 'x' || ch == 'X' || ch == '*') {
      str = read_digits(s, true, true);
      if (s.peek() == '.') {
        s >> ch;
        digit_str = read_digits(s, false, true);
      }
    } else
      s.putback(ch);
  } else {
    expn_before = false;
    s.putback(ch);
    str += read_digits(s, true, true);
    s >> ch;
    if (ch == '.') {
      digit_str = read_digits(s, false, true);
    } else 
      s.putback(ch);
  }

  if (!expn_before) {
    s >> ch;
    if (ch == 'e' || ch == 'E') {
      s >> std::ws;
      expn_str = read_digits(s);
    } else
      s.putback(ch);
  }

  // consume the ending comma or period, if any
  s >> ch;
  if (s.eof())
    s.clear();
  else if (ch != ';' && ch != ',') 
    s.putback(ch);

  construct(expn_str, str, digit_str);

  s.flags(old_flags);
  return !s.fail();
}

bool mp_real::read(const string &s) {
  std::istringstream is(s);
  return read(is);
}
