/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include <stdio.h>

#include <alloc.h>
#include <vector.h>
#include <hashtab.h>

#include "symfcn.h"


int part_itr_sz(vector *part)
{
  int i, e, tot;

  i = v_length(part) - 1;
  while (i >= 0 && v_elem(part, i) == 1)
    i--;
  
  if (i < 0)
    return 0;
  
  e = v_elem(part, i);
  tot = e + v_length(part) - i - 1;
  
  while (tot >= e-1)
    {
      v_elem(part, i++) = e-1;
      tot -= e-1;
    }
  if (tot > 0)
    v_elem(part, i++) = tot;
  v_length(part) = i;
  
  return 1;
}

int part_itr_sub(vector *part, vector *outer)
{
  int i, e;

  i = v_length(part) - 1;
  while (i >= 0 && v_elem(part, i) == 0)
    i--;
  if (i < 0)
    return 0;
  
  e = v_elem(part, i) - 1;
  if (e == 0)
    {
      v_length(part) = i;
      return 1;
    }
  
  while (i < v_length(outer))
    {
      int x = v_elem(outer, i);
      v_elem(part, i) = (e < x) ? e : x;
      i++;
    }
  v_length(part) = v_length(outer);
  
  return 1;
}

int part_itr_between(vector *part, vector *inner, vector *outer)
{
  int i, e;
  
  i = v_length(part) - 1;
  while (i >= 0 && v_elem(part, i) == v_elem(inner, i))
    i--;
  if (i < 0)
    return 0;
  
  e = v_elem(part, i) - 1;
  if (e == 0)
    {
      v_length(part) = i;
      return 1;
    }
  
  while (i < v_length(outer))
    {
      int x = v_elem(outer, i);
      v_elem(part, i) = (e < x) ? e : x;
      i++;
    }
  v_length(part) = v_length(outer);
  
  return 1;
}

int part_length(vector *p)
{
  int n = v_length(p);
  while (n > 0 && v_elem(p, n-1) == 0)
    n--;
  return n;
}

vector *part_conjugate(vector *p)
{
  int np, nc, j;
  vector *conj;
  
  np = v_length(p);
  if (np == 0)
    return v_new(0);
  
  nc = v_elem(p, 0);
  if (nc == 0)
    return v_new(0);
  
  conj = v_new(nc);
  j = 0;
  while (np > 0)
    {
      int jlim = v_elem(p, np - 1);
      while (j < jlim)
	v_elem(conj, j++) = np;
      np--;
    }
  
  return conj;
}

int part_subset(vector *p1, vector *p2)
{
  int n;
  
  n = part_length(p1);
  if (n > v_length(p2))
    return 0;
  
  for (n--; n >= 0; n--)
    if (v_elem(p1, n) > v_elem(p2, n))
      return 0;
  
  return 1;
}


void _chop_cols(vector *in1, vector *out)
{
  int n, m, mm, i;
  
  n = v_length(in1);
  m = v_elem(in1, n - 1);
  mm = v_elem(out, n - 1);
  if (m > mm)
    m = mm;
  
  for (i = 0; i < n; i++)
    {
      v_elem(in1, i) -= m;
      v_elem(out, i) -= m;
    }
  
  v_length(in1) = part_length(in1);
  v_length(out) = part_length(out);
}

void _chop_rows(vector *in1, vector *out)
{
  int m, i;
  
  m = 1;
  while (m < v_length(in1) && v_elem(in1, m) == v_elem(out, m))
    m++;
  
  for (i = 0; i < v_length(in1) - m; i++)
    v_elem(in1, i) = v_elem(in1, i + m);
  v_length(in1) -= m;
  
  for (i = 0; i < v_length(out) - m; i++)
    v_elem(out, i) = v_elem(out, i + m);
  v_length(out) -= m;
}


long long lrcoef(vector *outer, vector *inner1, vector *inner2)
{
  vector *out, *in1, *in2, *out_conj, *cont, *vtmp;
  int w_out, w_in1, w_in2, do_swap;
  int rows, i, j, stack_sz, sp, x;
  int *skewtab, *itab, *jtab, *max_tab;
  long long res;
  
  out = v_new_copy(outer);
  in1 = v_new_copy(inner1);
  in2 = v_new_copy(inner2);
  
  v_length(out) = part_length(out);
  v_length(in1) = part_length(in1);
  v_length(in2) = part_length(in2);
  
  while (v_length(out) > 0 &&
	 v_length(out) >= v_length(in1) &&
	 v_length(out) >= v_length(in2))
    {
      if (v_length(in1) == v_length(out))
	{
	  _chop_cols(in1, out);
	}
      else if (v_length(in2) == v_length(out))
	{
	  _chop_cols(in2, out);
	}
      else if (v_length(in1) > 0 && 
	       v_elem(in1, 0) == v_elem(out, 0))
	{
	  _chop_rows(in1, out);
	}
      else if (v_length(in2) > 0 && 
	       v_elem(in2, 0) == v_elem(out, 0))
	{
	  _chop_rows(in2, out);
	}
      else
	break;
    }
  
  w_in1 = v_sum(in1);
  w_in2 = v_sum(in2);
  w_out = v_sum(out);
  
  if ((! part_subset(in1, out)) || 
      (! part_subset(in2, out)) ||
      (w_out != w_in1 + w_in2))
    {
      v_free(in1);
      v_free(in2);
      v_free(out);
      return 0;
    }
  
  if (w_in1 <= 1 || w_in2 <= 1)
    {
      v_free(in1);
      v_free(in2);
      v_free(out);
      return 1;
    }
  
  out_conj = part_conjugate(out);
  
  do_swap = (w_in1 < w_in2);
  
  if (w_in1 == w_in2)
    {
      int c1 = v_elem(in1, 0);
      int r1 = v_length(in1);
      int n1 = (c1 < r1) ? c1 : r1;
      
      int c2 = v_elem(in2, 0);
      int r2 = v_length(in2);
      int n2 = (c2 < r2) ? c2 : r2;
      
      do_swap = (n1 < n2);
    }
  
  if (do_swap)
    {
      int tmp;
      
      tmp = w_in1;
      w_in1 = w_in2;
      w_in2 = tmp;
      
      vtmp = in1;
      in1 = in2;
      in2 = vtmp;
    }
  
  if (v_length(in2) > v_elem(in2, 0))
    {
      vtmp = out;
      out = out_conj;
      out_conj = vtmp;
      
      vtmp = part_conjugate(in1);
      v_free(in1);
      in1 = vtmp;
      
      vtmp = part_conjugate(in2);
      v_free(in2);
      in2 = vtmp;
    }
  
  rows = v_length(out);
  
  vtmp = v_new(rows);
  for (i = 0; i < v_length(in1); i++)
    v_elem(vtmp, i) = v_elem(in1, i);
  for ( ; i < rows; i++)
    v_elem(vtmp, i) = 0;
  v_free(in1);
  in1 = vtmp;
  
  cont = v_new_zero(v_length(in2));
  
  stack_sz = w_out - w_in1;
  skewtab = amalloc(stack_sz * sizeof(int));
  max_tab = amalloc(stack_sz * sizeof(int));
  
  itab = amalloc(stack_sz * sizeof(int));
  jtab = amalloc(stack_sz * sizeof(int));
  sp = 0;
  for (i = 0; i < v_length(out); i++)
    {
      int left = (i < v_length(in1)) ? v_elem(in1, i) : 0;
      for (j = v_elem(out, i) - 1; j >= left; j--)
        {
	  itab[sp] = i;
	  jtab[sp] = j;
	  sp++;
	}
    }
  
  res = 0;
  
  skewtab[0] = 0;
  v_elem(cont, 0) = 1;
  
  sp = 1;
  i = itab[1];
  j = jtab[1];
  max_tab[1] = (i == itab[0]) ? 0 : 
    v_length(in2) + i - v_elem(out_conj, j);;
  x = (j == jtab[0]) ? 1 : 0;
  
  while (sp > 0)
    {
      /* in each loop, do one of the following:
       *
       *   a) incr x
       *
       *   b) put x; incr sp; get x;
       *
       *   c) decr sp; get x; incr x;
       *
       * where:  put/get involves updating cont and skewtab and
       *
       *         incr sp includes setting max_tab[sp].
       */
      
      if (x > max_tab[sp])
	{
	  sp--;
	  
	  x = skewtab[sp];
	  v_elem(cont, x)--;
	  
	  x++;
	  continue;
	}
      
      if (v_elem(cont, x) == v_elem(in2, x) ||
	  (x > 0 && v_elem(cont, x) >= v_elem(cont, x-1)))
	{
	  x++;
	  continue;
	}
      
      if (sp + 1 == stack_sz)
	{
	  /* special case: incr res; then same as c) */
	  res++;
	  
	  sp--;
	  
	  x = skewtab[sp];
	  v_elem(cont, x)--;
	  
	  x++;
	  continue;
	}
      
      skewtab[sp] = x;
      v_elem(cont, x)++;
      
      sp++;
      i = itab[sp];
      j = jtab[sp];
      
      if (j < v_elem(out, i) - 1)
	max_tab[sp] = skewtab[sp - 1];
      else
	max_tab[sp] = v_length(in2) + i - v_elem(out_conj, j);
      
      x = 0;
      if (i > 0 && j >= v_elem(in1, i-1))
	x = skewtab[sp + v_elem(in1, i-1) - v_elem(out, i)] + 1;
    }
  
  afree(skewtab);
  afree(max_tab);
  afree(itab);
  afree(jtab);

  v_free(in1);
  v_free(in2);
  v_free(out);
  v_free(out_conj);
  v_free(cont);
  
  return res;
}


void st_setmin(skewtab *st, int y, int x)
{
  while (y < st->rows)
    {
      while (x >= v_elem(st->inner, y))
	{
	  int e;
	  if (y == 0 || x < v_elem(st->inner, y-1))
	    e = 0;
	  else
	    e = st->matrix[x + (y-1) * st->cols] + 1;

	  st->matrix[x + y * st->cols] = e;

	  v_elem(st->conts, e)++;

	  x--;
	}
      y++;

      if (y < st->rows)
	x = v_elem(st->outer, y) - 1;
    }
}


skewtab *st_new(vector *outer, vector *inner, vector *conts, int maxrows)
{
  int rows, cols, clen;
  skewtab *st;

  rows = v_length(outer);
  cols = (rows == 0) ? 0 : v_elem(outer, 0);

  st = amalloc(sizeof(skewtab) + sizeof(int) * (rows * cols - 1));
  st->outer = outer;
  st->inner = inner;
  
  clen = (conts != NULL) ? v_length(conts) : 0;
  clen += rows;
  st->conts = v_new_zero(clen);
  if (conts != NULL)
    {
      int i;
      for (i = 0; i < v_length(conts); i++)
	v_elem(st->conts, i) = v_elem(conts, i);
    }
  
  st->rows = rows;
  st->cols = cols;
  
  st_setmin(st, 0, v_elem(outer, 0) - 1);
  
  st->conjugate = NULL;
  if (maxrows >= clen)
    maxrows = 0;
  st->maxrows = maxrows;
  
  if (maxrows == 0)
    return st;
  
  if (v_elem(st->conts, maxrows) != 0)
    {
      st_free(st);
      return NULL;
    }
  
  st->conjugate = part_conjugate(outer);
  
  return st;
}

int st_next(skewtab *st)
{
  int x, y, e;

  for (y = st->rows - 1; y >= 0; y--)
    {
      int xlim = v_elem(st->outer, y);

      for (x = v_elem(st->inner, y); x < xlim; x++)
	{
	  int elim = (st->maxrows == 0)
	    ? v_length(st->conts) - 1
	    : st->maxrows + y - v_elem(st->conjugate, x);
	  
	  if (x != xlim - 1)
	    {
	      int next_cell = st->matrix[x+1 + y * st->cols];
	      if (next_cell < elim)
		elim = next_cell;
	    }
	  
	  e = st->matrix[x + y * st->cols];
	  v_elem(st->conts, e)--;

	  e++;
	  while (e <= elim &&
		 v_elem(st->conts, e) == v_elem(st->conts, e-1))
	    e++;

	  if (e <= elim)
	    {
	      st->matrix[x + y * st->cols] = e;
	      v_elem(st->conts, e)++;
	      st_setmin(st, y, x-1);
	      return 1;
	    }
	}
    }

  return 0;
}

void st_print(skewtab *st)
{
  int x, y;
  
  vector *inner = st->inner;
  vector *outer = st->outer;
  
  for (y = 0; y < v_length(outer); y++)
    {
      for (x = 0; x < v_elem(outer, y); x++)
	{
	  if (x < v_elem(inner, y))
	    putchar(' ');
	  else
	    printf("%d", st->matrix[x + y * st->cols]);
	}
      putchar('\n');
    }
}

void st_free(skewtab *st)
{
  v_free(st->inner);
  v_free(st->outer);
  v_free(st->conts);
  if (st->conjugate != NULL)
    v_free(st->conjugate);
  afree(st);
}


hashtab *skew(vector *outer, vector *inner, int maxrows)
{
  vector *out0, *in0;
  hashtab *res;
  vector *part;
  skewtab *st;
  int n, i;
  
  res = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  
  n = v_length(outer);
  if (v_length(inner) > n)
    return res;
  
  out0 = v_new_copy(outer);
  
  in0 = v_new_zero(n);
  for (i = 0; i < v_length(inner); i++)
    v_elem(in0, i) = v_elem(inner, i);
  
  if (! v_lesseq(in0, out0))
    {
      v_free(in0);
      v_free(out0);
      return res;
    }
  
  st = st_new(out0, in0, NULL, maxrows);
  
  part = v_new(n);
  
  do {
    int i, *valuep;
    
    for (i = 0; i < v_length(st->conts) && v_elem(st->conts, i) != 0; i++)
      v_elem(part, i) = v_elem(st->conts, i);
    v_length(part) = i;
    
    valuep = hash_mkfindint(res, part);
    if (hash_key_used)
      part = v_new(n);
    (*valuep)++;
  } while (st_next(st));
  
  st_free(st);
  v_free(part);
  
  return res;
}


hashtab *mult(vector *sh1, vector *sh2, int maxrows)
{
  skewtab *st;
  vector *part, *out0, *in0;
  hashtab *res;
  int n;
  
  res = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  
  if (v_sum(sh1) > v_sum(sh2))
    {
      vector *tmp = sh1;
      sh1 = sh2;
      sh2 = tmp;
    }
  
  out0 = v_new_copy(sh1);
  in0 = v_new_zero(v_length(sh1));
  st = st_new(out0, in0, sh2, maxrows);
  if (st == NULL)
    return res;
  
  n = v_length(sh1) + v_length(sh2);
  part = v_new(n);
  
  do {
    int i, *valuep;
    
    for (i = 0; i < v_length(st->conts) && v_elem(st->conts, i) != 0; i++)
      v_elem(part, i) = v_elem(st->conts, i);
    v_length(part) = i;
    
    valuep = hash_mkfindint(res, part);
    if (hash_key_used)
      part = v_new(n);
    (*valuep)++;
    
  } while (st_next(st));
  
  st_free(st);
  v_free(part);
  
  return res;
}


hashtab *schur_lc_mult(hashtab *lc1, hashtab *lc2, int maxrows)
{
  hash_itr itr1, itr2, itr;
  hashtab *pairs = hash_new((cmp_t) vp_cmp, (hash_t) vp_hash);
  hashtab *res;
  
  for (hash_first(lc1, itr1); hash_good(itr1); hash_next(itr1))
    {
      vector *v1 = hash_key(itr1);
      int c1 = hash_intvalue(itr1);
      
      for (hash_first(lc2, itr2); hash_good(itr2); hash_next(itr2))
	{
	  vector *v2 = hash_key(itr2);
	  int c2 = hash_intvalue(itr2);
	  
	  vecpair *vp = vp_new_unordered(v_new_copy(v1), v_new_copy(v2));
	  int *valp = hash_mkfind(pairs, vp);
	  *((int *) valp) += c1 * c2;
	  if (! hash_key_used)
	    vp_free(vp);
	}
    }
  
  res = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  
  for (hash_first(pairs, itr); hash_good(itr); hash_next(itr))
    {
      vecpair *vp = hash_key(itr);
      int c = hash_intvalue(itr);
      hashtab *prd;
      
      prd = mult(vp_first(vp), vp_second(vp), maxrows);
      lincomb_add_multiple(res, c, prd, (freekey_t) v_free1, NULL);
      hash_free(prd);  /* lincomb_add_multiple has frees the vectors */
    }
  
  free_vp_lincomb(pairs);
  
  return res;
}


hashtab *coprod(vector *part, int all)
{
  hashtab *res = hash_new((cmp_t) vp_cmp, (hash_t) vp_hash);
  int df, in_wt, wt = v_sum(part);
  vecpair *vp;
  vector *in = v_new_copy(part);
  
  do {
    in_wt = v_sum(in);
    
    if (2 * in_wt >= wt)
      {
	hashtab *sk = skew(part, in, 0);
	hash_itr itr;
	
	for (hash_first(sk, itr); hash_good(itr); hash_next(itr))
	  {
	    vector *in2 = hash_key(itr);
	    int coef = hash_intvalue(itr);
	    
	    df = 1;
	    if (2 * in_wt == wt && (df = v_cmp(in, in2)) < 0)
	      {
		v_free(in2);
		continue;
	      }
	    
	    vp = vp_new(v_new_copy(in), in2);
	    hash_insertint(res, vp, coef);
	    
	    if (all && df != 0)
	      {
		vp = vp_new(v_new_copy(in2), v_new_copy(in));
		hash_insertint(res, vp, coef);
	      }
	  }
	
	hash_free(sk);
      }
  } while (part_itr_sub(in, part));
  v_free(in);
  
  return res;
}


int rim_hook(vector *lambda, int rows, int cols, int *qp)
{
  int i, j, len, sign, q, n;
  
  len = v_length(lambda);
  n = rows + cols;

  q = 0;
  for (i = 0; i < len; i++)
    {
      int a = v_elem(lambda, i) + rows - i - 1;
      q += a / n;
      a %= n;
      v_elem(lambda, i) = a - rows + 1;
    }

  /* bubble sort :-( */
  sign = (rows & 1) ? 0 : q;
  for (i = 1; i < len; i++)
    {
      int a = v_elem(lambda, i);
      for (j = i; j > 0 && a > v_elem(lambda, j-1); j--)
	{
	  v_elem(lambda, j) = v_elem(lambda, j-1);
	}
      if (j > 0 && a == v_elem(lambda, j-1))
	return 0;
      v_elem(lambda, j) = a;
      sign += i - j;
    }
  
  for (i = 0; i < len; i++)
    {
      v_elem(lambda, i) += i;
      if (v_elem(lambda, i) < 0)
	return 0;
    }
  
  while (len > 0 && v_elem(lambda, len - 1) == 0)
    len--;
  v_length(lambda) = len;
  *qp = q;
  return (sign & 1) ? -1 : 1;
}

list *_quantum_reduce(hashtab* s, int rows, int cols)
{
  int sign, q, *valuep;
  list *qlist;
  hash_itr itr;
  
  qlist = l_newsz(10);
  
  for (hash_first(s, itr); hash_good(itr); hash_next(itr))
    {
      vector *lambda = hash_key(itr);
      int coef = hash_intvalue(itr);
      hashtab *tab;
      
      sign = rim_hook(lambda, rows, cols, &q);
      
      if (sign == 0)
	{
	  v_free(lambda);
	  continue;
	}
      
      while (q >= l_length(qlist))
	l_append(qlist, hash_new((cmp_t) v_cmp, (hash_t) v_hash));
      
      tab = l_elem(qlist, q);
      valuep = hash_mkfindint(tab, lambda);
      *valuep += sign * coef;
      if (! hash_key_used)
	v_free(lambda);
    }
  
  return qlist;
}

list *quantum_reduce(hashtab* s, int rows, int cols)
{
  list *qlist = _quantum_reduce(s, rows, cols);
  hash_free(s);
  return qlist;
}

void fusion_reduce(hashtab *lc, int rows, int cols, int opt_zero)
{
  int i, j, k, lamj;
  list *qlist = _quantum_reduce(lc, rows, cols);
  if (l_length(qlist) == 0)
    {
      hash_reset(lc);
      return;
    }
  hash_copy(lc, l_elem(qlist, 0));
  hash_free(l_elem(qlist, 0));
  for (i = 1; i < l_length(qlist); i++)
    {
      hashtab *tab = l_elem(qlist, i);
      hash_itr itr;
      
      for (hash_first(tab, itr); hash_good(itr); hash_next(itr))
	{
	  vector *lambda = hash_key(itr);
	  vector *mu;
	  
	  if (hash_intvalue(itr) == 0 && opt_zero == 0)
	    continue;
	  
	  mu = v_new(rows);
	  for (j = 0; j < rows; j++)
	    {
	      lamj = (j < v_length(lambda)) ? v_elem(lambda,j) : 0;
	      k = (j + i) % rows;
              v_elem(mu,k) = lamj + ((i+j) / rows) * cols + i;
	    }
	  hash_insertint(lc, mu, hash_intvalue(itr));
	}
      
      free_vec_lincomb(tab);
    }
  l_free(qlist);
}

