/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
extern int optind;

#include "vectarg.h"
#include "vector.h"

vector *get_vect_arg(int ac, char **av)
{
  int n, i, x;
  int *tmp;
  vector *res;
  char *endptr;
  char ch;
  
  if (optind == ac)
    return NULL;
  
  if (optind == 0)
    {
      optind++;
    }
  else
    {  
      /* skip any "-" or "/" argument */
      ch = *(av[optind]);
      if ((ch == '-' || ch == '/') && *(av[optind] + 1) == '\0')
	optind++;
    }
  
  tmp = amalloc((ac - optind) * sizeof(int));
  n = 0;
  
  while (optind < ac)
    {
      x = strtol(av[optind], &endptr, 10);
      if (endptr == av[optind] || *endptr != '\0')
	break;
      
      tmp[n++] = x;
      optind++;
    }
  
  if (n == 0)
    {
      afree(tmp);
      return NULL;
    }
  
  res = v_new(n);
  for (i = 0; i < n; i++)
    v_elem(res, i) = tmp[i];
  
  afree(tmp);
  
  return res;
}
