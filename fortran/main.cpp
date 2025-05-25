#include "config.h"

#define f_main FC_FUNC_(f_main, F_MAIN)

extern "C" void f_main();

int main() {
  f_main();
  return 0;
}

