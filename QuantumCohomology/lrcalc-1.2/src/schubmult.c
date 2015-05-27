/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
extern char *optarg;

#include "vectarg.h"
#include "lincomb.h"
#include "schublib.h"
#include "maple.h"

#define PROGNAME "schubmult"

  
void print_usage()
{
  fprintf(stderr, "usage: " PROGNAME " [-m] [-r rank] perm1 - perm2\n");
  exit(1);
}

int main(int ac, char **av)
{
  hashtab *s;
  vector *w1, *w2;
  int opt_maple = 0;
  int rank = 0;
  int c;
  
  if (setjmp(lrcalc_panic_frame))
    {
      fprintf(stderr, "out of memory.\n");
      exit(1);
    }

  while ((c = getopt(ac, av, "mr:")) != EOF)
    switch (c)
      {
      case 'm':
	opt_maple = 1;
	break;
      case 'r':
	rank = atoi(optarg);
	if (rank < 0)
	  print_usage();
	break;
      default:
	print_usage();
      }
  
  w1 = get_vect_arg(ac, av);
  w2 = get_vect_arg(ac, av);
  if (w1 == NULL || w2 == NULL)
    print_usage();
  
  s = mult_schubert(w1, w2, rank);
  
  if (opt_maple)
    maple_print_lincomb(s, "X", 1);
  else
    print_vec_lincomb(s, 0);
  
  v_free(w1);
  v_free(w2);
  free_vec_lincomb(s);
  
  memory_report;
  
  return 0;
}
