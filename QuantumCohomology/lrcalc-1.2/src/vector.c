/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include <stdio.h>
#include <stdarg.h>

#include "vector.h"


vector *v_new(int length)
{
  vector *v = (vector *) amalloc(sizeof(vector) + sizeof(int) * (length - 1));
  v_length(v) = length;
  return v;
}

vector *v_new_init(int length, ...)
{
  vector *v;
  va_list ap;
  int i;
  
  v = v_new(length);
  va_start(ap, length);
  for (i = 0; i < length; i++)
    v_elem(v, i) = va_arg(ap, int);
  va_end(ap);
  
  return v;
}

vector *v_new_zero(int length)
{
  vector *v = v_new(length);
  v_set_zero(v);
  return v;
}

vector *v_new_copy(vector *v)
{
  vector *n = v_new(v_length(v));
  v_copy(n, v);
  return n;
}

void v_free1(vector *v)
{
  v_free(v);
}


int v_sum(vector *v)
{
  int i, res;
  res = 0;
  for (i = 0; i < v_length(v); i++)
    res += v_elem(v, i);
  return res;
}

int v_lesseq(vector *v1, vector *v2)
{
  int i;
  i = 0;
  while (i < v_length(v1) && v_elem(v1, i) <= v_elem(v2, i))
    i++;
  return (i == v_length(v1));
}

int i_gcd(int x, int y)
{
  while (x != 0)
    {
      int z = y % x;
      y = x;
      x = z;
    }
  return abs(y);
}

int v_gcd(vector *v)
{
  size_t i;
  int res = 0;
  for (i = 0; i < v_length(v); i++)
    res = i_gcd(res, v_elem(v, i));
  return res;
}

int v_cmp(vector *v1, vector *v2)
{
  int i;

  if (v_length(v1) != v_length(v2))
    return v_length(v1) - v_length(v2);

  for (i = 0; i < v_length(v1); i++)
    if (v_elem(v1, i) != v_elem(v2, i))
      return v_elem(v1, i) - v_elem(v2, i);

  return 0;
}

void v_mult(vector *dst, int c, vector *src)
{
  size_t i;
  for (i = 0; i < v_length(dst); i++)
    v_elem(dst, i) = c * v_elem(src, i);
}

void v_div(vector *dst, vector *src, int c)
{
  size_t i;
  for (i = 0; i < v_length(dst); i++)
    v_elem(dst, i) = v_elem(src, i) / c;
}

int v_max(vector *v)
{
  int i, max, n;
  
  n = v_length(v);
  if (n == 0)
    return 0;
  
  max = v_elem(v, 0);
  for (i = 1; i < v_length(v); i++)
    if (max < v_elem(v, i))
      max = v_elem(v, i);
  
  return max;
}

hashkey_t v_hash(vector *v)
{
  int i;
  hashkey_t h, g;
  
  h = v_length(v);
  for (i = 0; i < v_length(v); i++)
    {
      h = (h << 3) + v_elem(v, i);
      g = h & 0xf0000000;
      h ^= (g >> 24);
      /* h ^= g; */
    }
  
  return h;
}


int v_elm_cmp(const void *x, const void *y)
{
  return *((int *)x) - *((int *)y);
}


void v_print(vector *v)
{
  int i;
  putchar('(');
  for (i = 0; i < v_length(v); i++)
    printf(i == 0 ? "%d" : ", %d", v_elem(v, i));
  putchar(')');
}

void v_printnl(vector *v)
{
  v_print(v);
  putchar('\n');
}


void print_vec_list(list *lst)
{
  size_t i;
  for (i = 0; i < l_length(lst); i++)
    v_printnl(l_elem(lst, i));
}

void free_vec_list(list *lst)
{
  int i;
  for (i = 0; i < l_length(lst); i++)
    v_free(l_elem(lst, i));
  l_free(lst);
}

void print_vec_set(set *s)
{
  vector *v;
  set_itr itr;
  for (v = s_first(s, itr); v != NULL; v = s_next(itr))
    v_printnl(v);
}

void free_vec_set(set *s)
{
  set_itr itr;
  vector *v;
  for (v = s_first(s, itr); v != NULL; v = s_next(itr))
    v_free(v);
  s_free(s);
}

void print_vec_lincomb(hashtab *ht, int opt_zero)
{
  hash_itr itr;
  for (hash_first(ht, itr); hash_good(itr); hash_next(itr))
    {
      if (hash_intvalue(itr) == 0 && opt_zero == 0)
	continue;
      printf("%d  ", hash_intvalue(itr));
      v_printnl(hash_key(itr));
    }
}

void free_vec_lincomb(hashtab *ht)
{
  hash_itr itr;
  for (hash_first(ht, itr); hash_good(itr); hash_next(itr))
    v_free(hash_key(itr));
  hash_free(ht);
}


vecpair *vp_new(vector *v1, vector *v2)
{
  vecpair *vp = amalloc(sizeof(vecpair));
  vp_first(vp) = v1;
  vp_second(vp) = v2;
  return vp;
}

vecpair *vp_new_unordered(vector *v1, vector *v2)
{
  if (v_cmp(v1, v2) <= 0)
    return vp_new(v1, v2);
  else
    return vp_new(v2, v1);
}

void vp_free(vecpair *vp)
{
  v_free(vp_first(vp));
  v_free(vp_second(vp));
  afree(vp);
}

int vp_cmp(vecpair *vp1, vecpair *vp2)
{
  int res = v_cmp(vp_first(vp1), vp_first(vp2));
  return res ? res : v_cmp(vp_second(vp1), vp_second(vp2));
}

hashkey_t vp_hash(vecpair *vp)
{
  hashkey_t h1 = v_hash(vp_first(vp));
  hashkey_t h2 = v_hash(vp_second(vp));

  return h2 ^ (h1 >> 16) ^ (h1 << 16);
}

void print_vp_lincomb(hashtab *ht)
{
  hash_itr itr;
  for (hash_first(ht, itr); hash_good(itr); hash_next(itr))
    {
      vecpair *vp = hash_key(itr);
      printf("%d  ", hash_intvalue(itr));
      v_print(vp_first(vp));
      printf("  ");
      v_printnl(vp_second(vp));
    }
}

void free_vp_lincomb(hashtab *ht)
{
  hash_itr itr;
  for (hash_first(ht, itr); hash_good(itr); hash_next(itr))
    vp_free(hash_key(itr));
  hash_free(ht);
}



list *find_extreme_vectors(list *veclist, int take_max)
{
  list *res = l_new();
  list *tmp = l_new_copy(veclist);
  
  while (v_length(tmp) > 0)
    {
      vector *x = l_deletelast(tmp);
      int i;
      
      for (i = 0; i < l_length(res); i++)
	{
	  int bad = take_max 
	    ? v_lesseq(x, l_elem(res, i))
	    : v_lesseq(l_elem(res, i), x);
	  
	  if (bad)
	    break;
	}
      
      if (i < l_length(res))
	continue;
      
      i = 0;
      while (i < l_length(tmp))
	{
	  int bad = take_max 
	    ? v_lesseq(x, l_elem(tmp, i))
	    : v_lesseq(l_elem(tmp, i), x);
	  
	  if (bad)
	    x = l_fastdelete(tmp, i);
	  else
	    i++;
	}
      
      l_append(res, x);
    }
  
  l_free(tmp);
  
  return res;
}
