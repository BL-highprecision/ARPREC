#include <stdlib.h>
#include <stdio.h>
#include <arprec/c_mp.h>

#define MALLOC(type, n)  ((type *) malloc((n) * sizeof(type)))

void print_mp(const char *str, const double *a) {
  if (str) printf("%s = ", str);
  c_mpout(a);
  printf("\n");
}

int main() {
  int nr_words, nr_digits;
  double *a, *b, *c;

  nr_digits = 291;
  nr_words = c_mpinit(nr_digits + 5);
  c_mpsetoutputprec(nr_digits);

  a = c_mpnew();
  b = c_mpnew();
  c = c_mpnew();

  c_mppi(a);

  print_mp("pi", a);

  c_mpeq_str("0.7", a);
  c_mpsin(a, b);
  c_mpcos(a, c);

  print_mp("sin(0.7)", b);
  print_mp("cos(0.7)", c);

  c_mpmul(b, b, a);
  c_mpmul(c, c, b);
  c_mpadd(a, b, c);
  c_mpsub_d(c, 1.0, a);

  print_mp("sin^2(0.7) + cos^2(0.7) - 1.0", a);

  c_mpfree(a);
  c_mpfree(b);
  c_mpfree(c);
  return 0;
}
