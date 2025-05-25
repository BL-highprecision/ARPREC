#ifndef PSLQ3_TEMPLATES_CC
#define PSLQ3_TEMPLATES_CC

template <class T>
int check_size_pslq(const matrix<T> &a, const matrix<T> &b, 
                    const T &mx1, const T &mx2) {
  int result = RESULT_CONTINUE;
  T u = 0.0;
  T t;

  matrix_min(a, u);
  if (u > mx1) {
    result = (u > mx2) ? RESULT_VERY_LARGE_VALUE : RESULT_LARGE_VALUE;
  }

  return result;
}

#endif
