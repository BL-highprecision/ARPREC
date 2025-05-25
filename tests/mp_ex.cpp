#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include <arprec/mp_int.h>

using std::cin;
using std::cout;
using std::endl;

int main()
{
  mp::mp_init(2005);

  mp_real a, b, s;


  // Character-to-MP assignment, generic MPREAL function, pi and e.
  // mp_real x = "1.234567890 1234567890 1234567890 e-100";
  // mp_real y = "1.234567890 1234567890 1234567890 e-100";
  double dx= 4.0, dy = 3., ds; // dy = 1./3.;
  mp_real x;
  mp_complex z;
  mp_int ia, ib;

//  cin >> x;
//  cout << "x = " << x << endl;

  mp::mpsetprecwords(10);
  mp_real y;

  cout << "input x:";
  cin  >> x;
  cout << "x = " << x << endl;

  printf(" dx = %22.18e\tdy = %22.18e\n", dx, dy);
  x = dx;
  y = dy;
  s = x * y;
  cout << "s = " << s;
  cout << "log(s) = " << log(s) << endl;
  cout << "-log(s) = " << -log(s) << endl;

  b = atan2(x, y);
  ds = fmod(dx, dy);
  b = -s;
  cout << "b = " << b << endl;
  cout << "ds = " << ds << endl;

  z = mp_complex(1.0, 5.0);
  cout << "z.real " << z.real << endl;
  cout << "z.imag " << z.imag << endl;
  z = x - z;
  cout << "z.real " << z.real << endl;
  cout << "z.imag " << z.imag << endl;

  
  ia = 3;
  ib = 4;
  ib = ia + ib;
  cout << "ia = " << ia << "ib = " << ib << endl;
  
  mp::mp_finalize();
  return 0;
}
