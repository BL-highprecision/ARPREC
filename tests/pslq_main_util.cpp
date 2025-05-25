#include <cstdio>
#include <arprec/mp_real.h>
#include "pslq1.h"

void print_usage() {
  puts("pslq1 [-d NUM] [-e NUM] [-h] [-m NUM] [-r NUM] [-s NUM] [-v]");

  puts("   -d NUM, --digits NUM");
  puts("           Use NUM digits of precision.");

  puts("   -e NUM, --eps NUM");
  puts("           Specify the base 10 log of the detection epsilon.");
  puts("           Usually set to  20 - ndp  or so, where ndp is the numebr of digits");
  puts("           of precision used.  Must no exceed the precision of input data.");

  puts("   -h, --help");
  puts("           Print out this help message.");

  puts("   -m NUM, --mode NUM");
  puts("           Select the data.  Mode can be any of the following:");
  puts("           0   Runs PSLQ on [ 1, a, a^2, ..., a^{n-1} ], where");
  puts("                 a = 3^(1/r) - 2^(1/s)  ");
  puts("               is an algebraic number of degree n = r * s.");

  puts("   -r NUM  Specify the integer r used in mode 0.");
  puts("   -s NUM  Specify the integer s used in mode 0.");

  puts("   -v, --verbose");
  puts("           Verbose mode.  More instances of -v will give increasing");
  puts("           verbosity.  Up to three levels are recognized.");
}

/* Clear all the timers. */
void clear_timers() {
  for (int i = 0; i < NR_TIMERS; i++)
    timers[i] = 0.0;
}

/* Report the time for the given timer. */
static void report(const char *name, int index) {
  printf("%20s   %10.2f         %6.2f%%\n", name, timers[index], 
         100.0 * timers[index] / timers[TIMER_PSLQ_TOTAL]);
}

/* Report all timers. */
void report_timers() {
  printf("                            Time (s)    Percentage\n");
  printf("--------------------------------------------------\n");
  report("MP Update:", TIMER_MP_UPDATE);
  report("MPM Update:", TIMER_MPM_UPDATE);
  report("MPM LQ Decomp:", TIMER_MPM_LQ);
  report("MP Init:", TIMER_MP_INIT);
  report("MPM Init:", TIMER_MPM_INIT);
  report("MPM Iterate:", TIMER_MPM_ITERATE);
  report("PSLQ Total:", TIMER_PSLQ_TOTAL);
}

/* Used by the PSLQ driver to construct the input to the PSLQ algorithm. */
void init_data(int mode, int n, int r, int s, 
               matrix<mp_real> &x, matrix<mp_real> &ans) {
  switch (mode) {
    case MODE_ALGEBRAIC_TEST: {
      if (n != r * s + 1) {
        fprintf(stderr, "ERROR: Total degree does not match.\n");
        exit(-1);
      }

      mp_real a, b;
      mp_real::mpnrt(mp_real(3.0), r, a, mp::prec_words);
      mp_real::mpnrt(mp_real(2.0), s, b, mp::prec_words);
      a -= b;
      b = a;
      x(0) = 1.0;
      for (int i = 1; i < n; i++, b *= a) {
        x(i) = b;
      }
      break;
    }

    default:
      fprintf(stderr, "ERROR: Invalid mode: %d.\n", mode);
      exit(-1);
  }
}

void parse_command(int argc, char **argv, int &mode, int &n, int &r, 
		               int &s, int &nr_digits, int &n_eps) {
	int val;
	char *arg;

  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-d") == 0 || strcmp(arg, "--digits") == 0) {
      if (++i < argc) {
        val = atoi(argv[i]);
        if (val < 40) {
          fprintf(stderr, "Number of digits must be larger than 40.  Using 40 digits ...\n");
          val = 40;
        }
        nr_digits = val;
      } else {
        fprintf(stderr, "A number must follow after -d or --digits.\n");
      }
    } else if (strcmp(arg, "-e") == 0 || strcmp(arg, "--eps") == 0) {
      if (++i < argc) {
        val = atoi(argv[i]);
        if (-val > nr_digits) {
          fprintf(stderr, "Detection epsilon must not exceed precision.\n");
        } else {
          n_eps = val;
        }
      } else {
        fprintf(stderr, "A number must follow after -e or --eps.\n");
      }
    } else if (strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0) {
      print_usage();
      exit(0);
    } else if (strcmp(arg, "-m") == 0 || strcmp(arg, "--mode") == 0) {
      if (++i < argc) {
        val = atoi(argv[i]);
        if (val < 0 || val >= NR_MODES) {
          fprintf(stderr, "Invalid mode: %d.  Using mode 0 ...\n", val);
        } else {
          mode = val;
        }
      } else {
        fprintf(stderr, "A number must follow after -m or --mode.\n");
      }
    } else if (strcmp(arg, "-r") == 0) {
      if (++i < argc) {
        val = atoi(argv[i]);
        if (val < 0) {
          fprintf(stderr, "r must be positive.  Using r = %d.\n", r);
        } else {
          r = val;
        }
      }
    } else if (strcmp(arg, "-s") == 0) {
      if (++i < argc) {
        val = atoi(argv[i]);
        if (val < 0) {
          fprintf(stderr, "s must be positive.  Using s = %d.\n", s);
        } else {
          s = val;
        }
      }
    } else if (strcmp(arg, "-v") == 0 || strcmp(arg, "--verbose") == 0) {
      debug_level++;
    } else {
      fprintf(stderr, "Ignoring unknown option: %s\n", arg);
    }
  }
}

