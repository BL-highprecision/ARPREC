#ifndef MP_MATRIX_H_
#define MP_MATRIX_H_

#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>

/* To enable index bounds check, set this flag to 1. */
#define MATRIX_DEBUG 0

/* 
 * class matrix<T>
 *
 * Simple template for a matrix class. 
 * Matrix elements are accessed through the (i, j) operator:  a(i, j).
 * All data are stored in column major format.  The access a(i) is
 * also allowed (for vectors).  If T is a type with contiguous data, 
 * then matrix<T> will store its elements in a contiguous block.
 *
 * Note that matrix<T> derives from matrix_base<T>, since 
 * The matrix_base<T> class does no initialization.  Use the matrix<T>
 * class instead.  This is done due to incompatible constructor 
 * between T = double and T = mp_real, to make the data contiguous.
 */

template <class T>
void print(const T &x);


template <class T> 
class matrix_base {
protected:
  int n, m;
  T *elements;

  matrix_base<T> () { }
  matrix_base<T> (matrix_base<T> &m) { }
public:

  inline matrix_base(int n0, int m0 = 1) : n(n0), m(m0) { }

  inline int size() const {
    return n * m;
  }

  inline void getSize(int &r, int &c) const {
    r = n;
    c = m;
  }

  inline T &elem(int i) const {
#if (MATRIX_DEBUG)
    if (i < 0 || i >= m*n) {
      std:cerr << "ERROR: matrix(i): index out of bounds: " << i << std::endl;
      exit(-1);
    } else
      return elements[i];
#else
    return elements[i];
#endif
  }

  inline T &elem(int i, int j) const {
#if (MATRIX_DEBUG)
    int k = n * j + i;
    if (k < 0 || k >= m*n) {
      std:cerr << "ERROR: matrix(i, j): index out of bounds: " 
        << i << ", " << j << std::endl;
      exit(-1);
    }
    return elements[n * j + i];
#else
    return elements[n * j + i];
#endif
  }

  inline T &operator()(int i) const {
    return elem(i);
  }

  inline T &operator()(int i, int j) const {
    return elem(i, j);
  }

  inline T *getElements() const {
    return elements;
  }

  void identity() {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        elem(i, j) = ((i == j) ? 1.0 : 0.0);
  }

  void zero() {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        elem(i, j) = 0.0;
  }

  void print(const char *name = NULL) const {
    if (name)
      std::cout << name << std::endl;
    if (m == 1) {
      for (int i = 0; i < n; i++) {
        std::cout << elem(i);
      }
      std::cout << std::endl;
    } else {
      for (int i = 0; i < n; i++) {
        std::cout << "row " << i << std::endl;
        for (int j = 0; j < n; j++) {
          std::cout << elem(i, j);
        }
        std::cout << std::endl;
      }
    }
  }

  inline matrix_base<T> &operator=(const matrix_base<T> &x) {
    if (x.m != m || x.n != n) {
      std::cerr << "ERROR (matrix<T>::operator=): dimension mismatch." << std::endl;
      exit(-1);
    }
    for (int i = 0; i < m*n; i++)
      elements[i] = x.elements[i];
    return *this;
  }

};

template <class T>
class matrix : public matrix_base<T> {
protected:
  matrix (matrix<T> &a) { }
public:
  matrix(int n0, int m0 = 1) : matrix_base<T>(n0, m0) {
    matrix_base<T>::elements = new T[m0 * n0];
  }

  ~matrix() {
    delete [] matrix_base<T>::elements;
  }
};

template <class T>
class mp_matrix : public matrix_base<T> {
public:
  mp_matrix(int n0, int m0 = 1) : matrix_base<T>(n0, m0) {
    int mn = matrix_base<T>::m * matrix_base<T>::n;
    matrix_base<T>::elements = T::alloc_array(mn);
  }

  ~mp_matrix() {
    T::free_array(matrix_base<T>::elements);
  }

  inline mp_matrix<T> &operator=(const mp_matrix<T> &x) {
    matrix_base<T>::operator=(x);
    return *this;
  }

};

template <>
class matrix<mp_real> : public mp_matrix<mp_real>  {
public:
  matrix(int n0, int m0 = 1) : mp_matrix<mp_real>(n0, m0) { }

  inline double *getData() const {
    return elements[0].mpr;
  }
};

template <>
class matrix<mp_complex> : public mp_matrix<mp_complex>  {
public:
  matrix(int n0, int m0 = 1) : mp_matrix<mp_complex>(n0, m0) { }

  inline double *getData() const {
    return elements[0].real.mpr;
  }
};

#endif

