/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include "vector.h"
#include "hashtab.h"
#include "claim.h"

#include "lincomb.h"
#include "schublib.h"

hashtab *trans(vector *w, int vars, hashtab *res)
{
  hashtab *tmp;
  hash_itr itr;
  vector *v;
  int n, nw, r, s, wr, ws, i, vi, vr, last;
  
  if (res == NULL)
    res = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  else
    hash_reset(res);
  
  nw = v_length(w);
  n = nw;
  while (n > 1 && v_elem(w, n-1) == n)
    n--;
  v_length(w) = n;
  
  r = n-1;
  while (r > 0 && v_elem(w, r-1) < v_elem(w, r))
    r--;
  if (r == 0)
    {
      vector *xx = v_new_zero(vars ? vars : 1);
      hash_insertint(res, xx, 1);
      v_length(w) = nw;
      return res;
    }
  if (vars < r)
    vars = r;
  
  s = r + 1;
  while (s < n && v_elem(w, r-1) > v_elem(w, s))
    s++;
  
  wr = v_elem(w, r-1);
  ws = v_elem(w, s-1);
  
  v = w;
  v_elem(v, s-1) = wr;
  v_elem(v, r-1) = ws;
  
  tmp = trans(v, vars, NULL);
  for (hash_first(tmp, itr); hash_good(itr); hash_next(itr))
    {
      vector *xx = hash_key(itr);
      v_elem(xx, r-1)++;
      hash_insertint(res, xx, hash_intvalue(itr));
    }
  
  last = 0;
  vr = v_elem(v, r-1);
  for (i = r-1; i >= 1; i--)
    {
      vi = v_elem(v, i-1);
      if (last < vi && vi < vr)
	{
	  last = vi;
	  v_elem(v, i-1) = vr;
	  v_elem(v, r-1) = vi;
	  trans(v, vars, tmp);
	  lincomb_add_multiple(res, 1, tmp, (freekey_t) v_free1, NULL);
	  v_elem(v, i-1) = vi;
	}
    }

  v_length(w) = nw;
  v_elem(w, s-1) = ws;
  v_elem(w, r-1) = wr;
  hash_free(tmp);
  
  return res;
}

#if 0
hashtab *monk(int i, hashtab *slc, int rank)
{
  hashtab *res = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  hash_itr itr;
  
  for (hash_first(slc, itr); hash_good(itr); hash_next(itr))
    {
      vector *u, *w = hash_key(itr);
      int j, t;
      int c = hash_intvalue(itr);
      int n = v_length(w);
      int wi = (i <= n) ? v_elem(w, i-1) : i;
      
      if (i <= n+1)
	{
	  int last, ulen = (i > n) ? i : n;
	  last = 0;
	  for (j = i-1; j >= 1; j--)
	    if (last < v_elem(w, j-1) && v_elem(w, j-1) < wi)
	      {
		last = v_elem(w, j-1);
		u = v_new(ulen);
		for (t = 0; t < ulen; t++)
		  v_elem(u, t) = v_elem(w, t);
		v_elem(u, j-1) = wi;
		v_elem(u, i-1) = last;
		
		if (lincomb_add_element(res, -c, u, (freekey_t) v_free1) == 0)
		  v_free(u);
	      }
	}
      else
	{
	  u = v_new(i);
	  for (t = 0; t < n; t++)
	    v_elem(u, t) = v_elem(w, t);
	  for (t = n; t < i-2; t++)
	    v_elem(u, t) = t+1;
	  v_elem(u, i-2) = i;
	  v_elem(u, i-1) = i-1;
	  
	  if (lincomb_add_element(res, -c, u, (freekey_t) v_free1) == 0)
	    v_free(u);
	}
      
      if (i >= n+1)
	{
	  u = v_new(i + 1);
	  for (t = 0; t < n; t++)
	    v_elem(u, t) = v_elem(w, t);
	  for (t = n; t < i; t++)
	    v_elem(u, t) = t+1;
	  v_elem(u, i-1) = i+1;
	  v_elem(u, i) = i;
	  
	  if (lincomb_add_element(res, c, u, (freekey_t) v_free1) == 0)
	    v_free(u);
	}
      else
	{
	  int last;
	  last = (((unsigned) -1) >> 1);
	  for (j = i+1; j <= n; j++)
	    if (wi < v_elem(w, j-1) && v_elem(w, j-1) < last)
	      {
		last = v_elem(w, j-1);
		u = v_new(n);
		for (t = 0; t < n; t++)
		  v_elem(u, t) = v_elem(w, t);
		v_elem(u, i-1) = last;
		v_elem(u, j-1) = wi;
		
		if (lincomb_add_element(res, c, u, (freekey_t) v_free1) == 0)
		  v_free(u);
	      }
	  if (last > n && n < rank)
	    {
	      u = v_new(n+1);
	      for (t = 0; t < n; t++)
		v_elem(u, t) = v_elem(w, t);
	      v_elem(u, i-1) = n+1;
	      v_elem(u, n) = wi;
	      
	      if (lincomb_add_element(res, c, u, (freekey_t) v_free1) == 0)
		v_free(u);
	    }
	}
    }
  
  return res;
}
#endif

void _monk_add(int i, hashtab *slc, int rank, hashtab *res);

hashtab *monk(int i, hashtab *slc, int rank)
{
  hashtab *res = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  _monk_add(i, slc, rank, res);
  return res;
}

void _monk_add(int i, hashtab *slc, int rank, hashtab *res)
{
  hash_itr itr;
  
  for (hash_first(slc, itr); hash_good(itr); hash_next(itr))
    {
      vector *u, *w = hash_key(itr);
      int j, t;
      int c = hash_intvalue(itr);
      int n = v_length(w);
      int wi = (i <= n) ? v_elem(w, i-1) : i;
      
      if (i <= n+1)
	{
	  int last, ulen = (i > n) ? i : n;
	  last = 0;
	  for (j = i-1; j >= 1; j--)
	    if (last < v_elem(w, j-1) && v_elem(w, j-1) < wi)
	      {
		last = v_elem(w, j-1);
		u = v_new(ulen);
		for (t = 0; t < ulen; t++)
		  v_elem(u, t) = v_elem(w, t);
		v_elem(u, j-1) = wi;
		v_elem(u, i-1) = last;
		
		if (lincomb_add_element(res, -c, u, (freekey_t) v_free1) == 0)
		  v_free(u);
	      }
	}
      else
	{
	  u = v_new(i);
	  for (t = 0; t < n; t++)
	    v_elem(u, t) = v_elem(w, t);
	  for (t = n; t < i-2; t++)
	    v_elem(u, t) = t+1;
	  v_elem(u, i-2) = i;
	  v_elem(u, i-1) = i-1;
	  
	  if (lincomb_add_element(res, -c, u, (freekey_t) v_free1) == 0)
	    v_free(u);
	}
      
      if (i >= n+1)
	{
	  u = v_new(i + 1);
	  for (t = 0; t < n; t++)
	    v_elem(u, t) = v_elem(w, t);
	  for (t = n; t < i; t++)
	    v_elem(u, t) = t+1;
	  v_elem(u, i-1) = i+1;
	  v_elem(u, i) = i;
	  
	  if (lincomb_add_element(res, c, u, (freekey_t) v_free1) == 0)
	    v_free(u);
	}
      else
	{
	  int last;
	  last = (((unsigned) -1) >> 1);
	  for (j = i+1; j <= n; j++)
	    if (wi < v_elem(w, j-1) && v_elem(w, j-1) < last)
	      {
		last = v_elem(w, j-1);
		u = v_new(n);
		for (t = 0; t < n; t++)
		  v_elem(u, t) = v_elem(w, t);
		v_elem(u, i-1) = last;
		v_elem(u, j-1) = wi;
		
		if (lincomb_add_element(res, c, u, (freekey_t) v_free1) == 0)
		  v_free(u);
	      }
	  if (last > n && n < rank)
	    {
	      u = v_new(n+1);
	      for (t = 0; t < n; t++)
		v_elem(u, t) = v_elem(w, t);
	      v_elem(u, i-1) = n+1;
	      v_elem(u, n) = wi;
	      
	      if (lincomb_add_element(res, c, u, (freekey_t) v_free1) == 0)
		v_free(u);
	    }
	}
    }
}


#if 0
hashtab *mult_poly_schubert(hashtab *poly, vector *perm, int rank)
{
  int maxvar, i, c;
  list *tmplist;
  hashtab *res, *res1, *poly1;
  hash_itr itr;
  
  if (hash_card(poly) == 0)
    return poly;
  
  maxvar = 0;
  for (hash_first(poly, itr); hash_good(itr); hash_next(itr))
    {
      vector *xx = hash_key(itr);
      c = hash_intvalue(itr);
      i = v_length(xx);
      if (c != 0 && i > maxvar)
	{
	  while (i > maxvar && v_elem(xx, i-1) == 0)
	    i--;
	  if (i > maxvar)
	    maxvar = i;
	}
    }
  
  if (maxvar == 0)
    {
      hash_first(poly, itr);
      c = hash_intvalue(itr);
      v_free(hash_key(itr));
      hash_reset(poly);   /* LOOK OUT FOR MEMORY LEAK!! */
      hash_insertint(poly, v_new_copy(perm), c);
      return poly;
    }
  
  tmplist = l_new();
  for (hash_first(poly, itr); hash_good(itr); hash_next(itr))
    {
      vector *xx = hash_key(itr);
      if (v_length(xx) >= maxvar && v_elem(xx, maxvar-1) != 0)
	{
	  l_append(tmplist, xx);
	  l_appendint(tmplist, hash_intvalue(itr));
	}
    }
  poly1 = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  for (i = 0; i < l_length(tmplist); i += 2)
    {
      vector *xx = l_elem(tmplist, i);
      c = l_intelem(tmplist, i+1);
      hash_remove(poly, xx);
      v_elem(xx, maxvar-1)--;
      while (v_length(xx) > 0 && v_elem(xx, v_length(xx)-1) == 0)
	v_length(xx)--;
      hash_insertint(poly1, xx, c);
    }
  l_free(tmplist);
  
  res1 = mult_poly_schubert(poly1, perm, rank);
  res = monk(maxvar, res1, rank);
  free_vec_lincomb(res1);
  
  if (hash_card(poly) != 0)
    {
      poly = mult_poly_schubert(poly, perm, rank);
      lincomb_add_multiple(res, 1, poly, (freekey_t) v_free1, NULL);
    }
  hash_free(poly);
  
  return res;
}
#endif


void _mult_ps(void **poly, int n, int maxvar, vector *perm, int rank,
	      hashtab *res);

hashtab *mult_poly_schubert(hashtab *poly, vector *perm, int rank)
{
  hash_itr itr;
  void **p;
  int i, j, n, maxvar, svlen;
  
  n = hash_card(poly);
  if (n == 0)
    return poly;
  
  p = (void **) amalloc(2 * n * sizeof(void *));
  i = 0;
  maxvar = 0;
  for (hash_first(poly, itr); hash_good(itr); hash_next(itr))
    {
      vector *xx = hash_key(itr);
      j = v_length(xx);
      while (j > 0 && v_elem(xx, j-1) == 0)
	j--;
      v_length(xx) = j;
      if (maxvar < j)
	maxvar = j;
      p[i++] = hash_key(itr);
      p[i++] = (void *) (long) hash_intvalue(itr);
    }
  claim(i == 2 * n);
  hash_reset(poly);
  
  svlen = v_length(perm);
  v_length(perm) = perm_group(perm);
  _mult_ps(p, n, maxvar, perm, rank, poly);
  v_length(perm) = svlen;
  
  for (i = 0; i < n; i++)
    v_free(p[2*i]);
  afree(p);
  
  return poly;
}

void _mult_ps(void **poly, int n, int maxvar, vector *perm, int rank,
	      hashtab *res)
{
  int i, j, c, lnxx, mv0, mv1;
  hashtab *res1;
  void *cc;
  
  if (maxvar == 0)
    {
      vector *w = v_new_copy(perm);  /* FIXME: OPTIMIZE! */
      c = (int) (long) poly[1];
      if (lincomb_add_element(res, c, w, (freekey_t) v_free1) == 0)
	v_free(w);
      return;
    }
  
  mv0 = 0;
  mv1 = 0;
  j = 0;
  for (i = 0; i < n; i++)
    {
      vector *xx = poly[2*i];
      lnxx = v_length(xx);
      if (lnxx < maxvar)
	{
	  if (mv0 < lnxx)
	    mv0 = lnxx;
	}
      else
	{
	  v_elem(xx, maxvar - 1)--;
	  while (lnxx > 0 && v_elem(xx, lnxx-1) == 0)
	    lnxx--;
	  v_length(xx) = lnxx;
	  if (mv1 < lnxx)
	    mv1 = lnxx;
	  poly[2*i] = poly[2*j];
	  poly[2*j] = xx;
	  cc = poly[2*i+1];
	  poly[2*i+1] = poly[2*j+1];
	  poly[2*j+1] = cc;
	  j++;
	}
    }
  
  res1 = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  _mult_ps(poly, j, mv1, perm, rank, res1);
  _monk_add(maxvar, res1, rank, res);
  free_vec_lincomb(res1);
  
  if (j < n)
    _mult_ps(poly + 2*j, n - j, mv0, perm, rank, res);
}


int num_inversions(vector *w)
{
  int i, j, n, res;
  n = v_length(w);
  res = 0;
  for (i = 0; i < n-1; i++)
    for (j = i+1; j < n; j++)
      if (v_elem(w, i) > v_elem(w, j))
	res++;
  return res;
}

int perm_group(vector *w)
{
  int i = v_length(w);
  while (i > 1 && v_elem(w, i-1) == i)
    i--;
  return i;
}

hashtab *mult_schubert(vector *w1, vector *w2, int rank)
{
  hashtab *poly;
  int svlen1, svlen2, w1len, w2len;
  
  svlen1 = v_length(w1);
  v_length(w1) = perm_group(w1);
  w1len = num_inversions(w1);
  
  svlen2 = v_length(w2);
  v_length(w2) = perm_group(w2);
  w2len = num_inversions(w2);
  
  if (rank == 0)
    {
      rank = (((unsigned) -1) >> 1);
    }
  else
    {
      /* FIXME: one can say exactly when the product is zero. */
      
      if (2 * (w1len + w2len) > rank * (rank - 1)
	  || v_length(w1) > rank
	  || v_length(w2) > rank)
	{
	  v_length(w1) = svlen1;
	  v_length(w2) = svlen2;
	  return hash_new((cmp_t) v_cmp, (hash_t) v_hash);
	}
    }
  
  if (w1len <= w2len)
    {
      poly = trans(w1, 0, NULL);
      poly = mult_poly_schubert(poly, w2, rank);
    }
  else
    {
      poly = trans(w2, 0, NULL);
      poly = mult_poly_schubert(poly, w1, rank);
    }
  
  v_length(w1) = svlen1;
  v_length(w2) = svlen2;
  
  return poly;
}

list *all_strings(vector *dimvec)
{
  int n, ld, i, j, k;
  vector *str, *cntvec;
  list *res;
  
  ld = v_length(dimvec);
  cntvec = v_new_zero(ld);
  n = v_elem(dimvec, ld - 1);
  
  str = v_new(n);
  j = 0;
  for (i = 0; i < ld; i++)
    {
      while (j < v_elem(dimvec, i))
	v_elem(str, j++) = i;
    }
  
  res = l_new();  /* FIXME: calculate length first */
  while (1)
    {
      l_append(res, v_new_copy(str));
      
      j = n-1;
      v_elem(cntvec, v_elem(str, j))++;
      while (j > 0 && v_elem(str, j-1) >= v_elem(str, j))
	{
	  v_elem(cntvec, v_elem(str, --j))++;
	}
      if (j == 0)
	break;
      
      i = v_elem(str, j-1);
      v_elem(cntvec, i++)++;
      while(v_elem(cntvec, i) == 0)
	i++;
      v_elem(str, j-1) = i;
      v_elem(cntvec, i)--;
      
      for (i = 0; i < ld; i++)
	{
	  for (k = 0; k < v_elem(cntvec, i); k++)
	    v_elem(str, j++) = i;
	  v_elem(cntvec, i) = 0;
	}
    }
  
  v_free(cntvec);
  v_free(str);
  return res;
}

list *all_perms(int n)
{
  list *res;
  vector *dimvec = v_new(n);
  int i;
  for (i = 0; i < n; i++)
    v_elem(dimvec, i) = i + 1;
  res = all_strings(dimvec);
  v_free(dimvec);
  return res;
}

vector *string2perm(vector *str)
{
  int N, n, i, j;
  vector *dimvec, *perm;
  
  n = v_length(str);
  
  N = 0;
  for (i = 0; i < n; i++)
    if (N < v_elem(str, i))
      N = v_elem(str, i);
  N++;
  
  dimvec = v_new_zero(N);
  for (i = 0; i < n; i++)
    v_elem(dimvec, v_elem(str, i))++;
  for (i = 1; i < N; i++)
    v_elem(dimvec, i) += v_elem(dimvec, i-1);
  
  perm = v_new(n);
  
  for (i = n-1; i >= 0; i--)
    {
      j = v_elem(str, i);
      v_elem(perm, --v_elem(dimvec, j)) = i;
    }
  
  v_free(dimvec);
  return perm;
}

vector *str2dimvec(vector *str)
{
  int i, n;
  vector *res;
  n = 0;
  for (i = 0; i < v_length(str); i++)
    if (n <= v_elem(str, i))
      n = v_elem(str, i) + 1;
  res = v_new_zero(n);
  for (i = 0; i < v_length(str); i++)
    v_elem(res, v_elem(str, i))++;
  for (i = 1; i < n; i++)
    v_elem(res, i) += v_elem(res, i-1);
  return res;
}

vector *perm2string(vector *perm, vector *dimvec)
{
  int n, i, j;
  vector *res;
  
  n = v_elem(dimvec, v_length(dimvec) - 1);
  res = v_new(n);
  j = 0;
  for (i = 0; i < v_length(dimvec); i++)
    while (j < v_elem(dimvec, i))
      {
	int wj = (j < v_length(perm)) ? v_elem(perm, j) : j+1;
	v_elem(res, wj - 1) = i;
	j++;
      }
  
  return res;
}

hashtab *mult_str_schubert(vector *str1, vector *str2)
{
  vector *dv1, *dv2, *w1, *w2;
  hashtab *s, *res;
  hash_itr itr;
  int i;
  
  if (v_length(str1) != v_length(str2))
    return NULL;
  dv1 = str2dimvec(str1);
  dv2 = str2dimvec(str2);
  
  if (v_cmp(dv1, dv2) != 0)
    {
      v_free(dv1);
      v_free(dv2);
      return NULL;
    }
  v_free(dv2);
  
  w1 = string2perm(str1);
  w2 = string2perm(str2);
  /* fixme */
  for (i = 0; i < v_length(w1); i++)
    {
      v_elem(w1, i)++;
      v_elem(w2, i)++;
    }
  s = mult_schubert(w1, w2, v_length(w1));
  v_free(w1);
  v_free(w2);
  
  res = hash_new((cmp_t) v_cmp, (hash_t) v_hash);
  for (hash_first(s, itr); hash_good(itr); hash_next(itr))
    {
      vector *str = perm2string(hash_key(itr), dv1);
      hash_insertint(res, str, hash_intvalue(itr));
    }
  
  free_vec_lincomb(s);
  v_free(dv1);
  
  return res;
}

