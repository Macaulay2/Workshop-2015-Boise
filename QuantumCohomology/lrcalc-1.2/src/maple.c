/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include <stdio.h>
#include <vector.h>
#include <hashtab.h>
#include "maple.h"


void maple_print_term(int c, vector *v, char *letter)
{
  int x, i;
  
  putchar((c < 0) ? '-' : '+');
  c = abs(c);
  printf("%d*%s[", c, letter);
  
  for (i = 0; i < v_length(v); i++)
    {
      if (i > 0)
	putchar(',');
      x = v_elem(v, i);
      printf("%d", x);
    }
  putchar(']');
}

void maple_print_lincomb(hashtab *ht, char *letter, int nl)
{
  hash_itr itr;
  putchar('0');
  for (hash_first(ht, itr); hash_good(itr); hash_next(itr))
    {
      if (hash_intvalue(itr) == 0)
	continue;
      maple_print_term(hash_intvalue(itr), hash_key(itr), letter);
    }
  if (nl)
    putchar('\n');
}
