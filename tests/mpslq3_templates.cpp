#ifndef MPSLQ_TEMPLATES_CC
#define MPSLQ_TEMPLATES_CC

template <class T>
int update_mpslq2(const matrix<T> &da, const matrix<T> &db, 
		  matrix<mp_real> &a, matrix<mp_real> &b, 
                  matrix<mp_real> &h, matrix<mp_real> &y, 
		  const mp_real &eps, const mp_real &teps) {
    int result = update_mpslq(da, db, b, h, y, eps, teps);
    matmul_left(da, a);
    return result;
}

#endif
