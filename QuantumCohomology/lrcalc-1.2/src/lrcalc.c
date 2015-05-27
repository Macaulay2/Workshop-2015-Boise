/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
extern char *optarg;

#include <vectarg.h>

#include "symfcn.h"
#include "maple.h"

#define MULT_USAGE \
"lrcalc mult [-mz] [-r rows] [-q rows,cols] [-f rows,level] part1 - part2\n"
#define SKEW_USAGE \
"lrcalc skew [-m] [-r rows] outer / inner\n"
#define COPROD_USAGE \
"lrcalc coprod [-a] part\n"
#define LRCOEF_USAGE \
"lrcalc lrcoef outer - inner1 - inner2\n"
#define LRTAB_USAGE \
"lrcalc lrtab [-r rows] outer / inner\n"


void print_usage()
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, MULT_USAGE);
  fprintf(stderr, SKEW_USAGE);
  fprintf(stderr, COPROD_USAGE);
  fprintf(stderr, LRCOEF_USAGE);
  fprintf(stderr, LRTAB_USAGE);
  exit(1);
}


void mult_usage()
{
  fprintf(stderr, "Usage: " MULT_USAGE);
  exit(1);
}

int mult_main(int ac, char **av)
{
  hashtab *s;
  vector *sh1, *sh2;
  int c, wt1, wt2;
  int opt_maple = 0;
  int opt_zero = 0;
  int opt_rows = 0;
  int opt_cols = 0;
  int opt_quantum = 0;
  int opt_fusion = 0;
  char *p;
  
  while ((c = getopt(ac, av, "mzr:q:f:")) != EOF)
    switch (c)
      {
      case 'm':
	opt_maple = 1;
	break;
      case 'z':
	opt_zero = 1;
	break;
      case 'r':
	opt_rows = atoi(optarg);
	if (opt_rows <= 0)
	  mult_usage();
	break;
      case 'q':
      case 'f':
	if (c == 'q')
	  opt_quantum = 1;
	else
	  opt_fusion = 1;
	opt_rows = strtol(optarg, &p, 10);
	if (p == NULL || *p != ',')
	  mult_usage();
	opt_cols = atoi(p + 1);
	if (opt_rows <= 0 || opt_cols <= 0)
	  mult_usage();
	break;
      default:
	mult_usage();
      }
  
  sh1 = get_vect_arg(ac, av);
  sh2 = get_vect_arg(ac, av);

  if (sh1 == NULL || sh2 == NULL)
    mult_usage();

  wt1 = v_sum(sh1);
  wt2 = v_sum(sh2);
  
  s = mult(sh1, sh2, opt_rows);
  
  if (opt_maple)
    printf("0");
  
  if (opt_quantum)
    {
      int n = opt_rows + opt_cols;
      list *qlist = quantum_reduce(s, opt_rows, opt_cols);  
      int i;
      
      for (i = 0; i < l_length(qlist); i++)
	{
	  hashtab *tab = l_elem(qlist, i);
	  char symbol[15];
	  sprintf(symbol, "q^%d*s", i);
	  
	  if (! opt_maple)
	    {
	      print_vec_lincomb(tab, opt_zero);
	    }
	  else if (wt1 + wt2 != n * i)
	    {
	      maple_print_lincomb(tab, symbol, 0);
	    }
	  else if (hash_card(tab) > 0)
	    {
	      hash_itr itr;
	      hash_first(tab, itr);
	      if (hash_intvalue(itr) != 0 || opt_zero)
		printf("%+d*q^%d", hash_intvalue(itr), i);
	    }
	  free_vec_lincomb(tab);
	}
      
      l_free(qlist);
      if (opt_maple)
	putchar('\n');
    }
  else
    {
      if (opt_fusion)
	fusion_reduce(s, opt_rows, opt_cols, opt_zero);
      
      if (opt_maple)
	maple_print_lincomb(s, "s", 1);
      else
	print_vec_lincomb(s, opt_zero);
      free_vec_lincomb(s);
    }

  v_free(sh1);
  v_free(sh2);
  
  memory_report;
  
  return 0;
}


void skew_usage()
{
  fprintf(stderr, "Usage: " SKEW_USAGE);
  exit(1);
}

int skew_main(int ac, char **av)
{
  hashtab *s;
  vector *outer, *inner;
  int c;
  int opt_maple = 0;
  int opt_rows = 0;

  while ((c = getopt(ac, av, "mr:")) != EOF)
    switch (c)
      {
      case 'm':
	opt_maple = 1;
	break;
      case 'r':
	opt_rows = atoi(optarg);
	break;
      default:
	skew_usage();
      }
  
  outer = get_vect_arg(ac, av);
  inner = get_vect_arg(ac, av);
  
  if (inner == NULL)
    skew_usage();
  
  s = skew(outer, inner, opt_rows);
  if (opt_maple)
    maple_print_lincomb(s, "s", 1);
  else
    print_vec_lincomb(s, 0);
  
  v_free(outer);
  v_free(inner);
  free_vec_lincomb(s);
  
  memory_report;
  
  return 0;
}


void coprod_usage()
{
  fprintf(stderr, "Usage: " COPROD_USAGE);
  exit(1);
}

int coprod_main(int ac, char **av)
{
  hashtab *s;
  vector *part;
  int c, opt_all = 0;
  
  if (setjmp(lrcalc_panic_frame))
    {
      fprintf(stderr, "out of memory.\n");
      exit(1);
    }
  
  while ((c = getopt(ac, av, "a")) != EOF)
    switch (c)
      {
      case 'a':
	opt_all = 1;
	break;
      default:
	coprod_usage();
      }
  
  part = get_vect_arg(ac, av);
  if (part == NULL)
    coprod_usage();
  
  s = coprod(part, opt_all);
  print_vp_lincomb(s);

  free_vp_lincomb(s);
  v_free(part);
  
  memory_report;
  
  return 0;
}


void lrcoef_usage()
{
  fprintf(stderr, "Usage: " LRCOEF_USAGE);
  exit(1);
}

int lrcoef_main(int ac, char **av)
{
  vector *lm, *mu, *nu;

  nu = get_vect_arg(ac, av);
  lm = get_vect_arg(ac, av);
  mu = get_vect_arg(ac, av);
  
  if (mu == NULL)
    lrcoef_usage();
  
  printf("%lld\n", lrcoef(nu, lm, mu));
  
  v_free(lm);
  v_free(mu);
  v_free(nu);
  
  memory_report;
  
  return 0;
}


void lrtab_usage()
{
  fprintf(stderr, "Usage: " LRTAB_USAGE);
  exit(1);
}

void print_lrskew_set(vector *outer, vector *inner, int opt_rows)
{
  vector *out0, *in0;
  skewtab *st;
  int n, i;

  n = v_length(outer);
  if (v_length(inner) > n)
    return;
  
  out0 = v_new_copy(outer);
  
  in0 = v_new_zero(n);
  for (i = 0; i < v_length(inner); i++)
    v_elem(in0, i) = v_elem(inner, i);
  
  if (! v_lesseq(in0, out0))
    {
      v_free(in0);
      v_free(out0);
    }
  
  st = st_new(out0, in0, NULL, opt_rows);
  do {
    st_print(st);
    putchar('\n');
  } while (st_next(st));

  st_free(st);
}

int lrtab_main(int ac, char **av)
{
  int c;
  vector *outer, *inner;
  int opt_rows = 0;

  if (setjmp(lrcalc_panic_frame))
    {
      fprintf(stderr, "out of memory.\n");
      exit(1);
    }
  
  while ((c = getopt(ac, av, "mr:")) != EOF)
    switch (c)
      {
      case 'r':
	opt_rows = atoi(optarg);
	break;
      default:
	lrtab_usage();
      }
  
  outer = get_vect_arg(ac, av);
  inner = get_vect_arg(ac, av);
  
  if (inner == NULL)
    lrtab_usage();
  
  print_lrskew_set(outer, inner, opt_rows);
  
  v_free(outer);
  v_free(inner);
  
  memory_report;
  
  return 0;
}


int main(int ac, char **av)
{
  char *cmd;
  extern char ** environ;

  if (setjmp(lrcalc_panic_frame))
    {
      fprintf(stderr, "out of memory.\n");
      exit(1);
    }

  if (ac < 2)
    print_usage();

  cmd = av[1];
  if (cmd[0] == 'l' && cmd[1] == 'r')
    cmd += 2;

  if (strcmp(cmd, "mult") == 0)
    mult_main(ac-1, av+1);
  else if (strcmp(cmd, "skew") == 0)
    skew_main(ac-1, av+1);
  else if (strcmp(cmd, "coprod") == 0)
    coprod_main(ac-1, av+1);
  else if (strcmp(cmd, "coef") == 0)
    lrcoef_main(ac-1, av+1);
  else if (strcmp(cmd, "tab") == 0)
    lrtab_main(ac-1, av+1);
  else
    print_usage();

  return 0;
}

