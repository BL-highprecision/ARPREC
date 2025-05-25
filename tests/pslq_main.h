#ifndef PSLQ_MAIN_H
#define PSLQ_MAIN_H

void print_usage();
void parse_command(int argc, char **argv, int &mode, int &n, int &r, 
		               int &s, int &nr_digits, int &n_eps);

#endif
