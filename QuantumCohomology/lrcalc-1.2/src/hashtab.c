/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include "hashtab.h"
#include "vector.h"

int hash_key_used;
void *hash_removed_key;

hashtab *hash_new(cmp_t cmp, hash_t hsh)
{
  hashtab *s = amalloc(sizeof(hashtab));
  hash_cmp(s) = cmp;
  hash_hash(s) = hsh;

  hash_tabsz(s) = INIT_TABSZ;
  hash_table(s) = amalloc(INIT_TABSZ * sizeof(size_t));

  s->elts_sz = INIT_ELTSZ;
  hash_elts(s) = amalloc(INIT_ELTSZ * sizeof(struct _hashelt));

  hash_reset(s);

  return s;
}

void hash_free(hashtab *s)
{
  afree(hash_table(s));
  afree(hash_elts(s));
  afree(s);
}

void hash_reset(hashtab *s)
{
  size_t i;
  hash_card(s) = 0;
  for (i = 0; i < hash_tabsz(s); i++)
    hash_table(s)[i] = _S_END;
  for (i = 0; i < hash_eltsz(s); i++)
    s->elts[i].next = i + 1;
  s->elts[hash_eltsz(s) - 1].next = _S_END;
  s->free_elts = 0;
}

void hash_makeroom(hashtab *s, size_t sz);

void hash_copy(hashtab *dst, hashtab *src)
{
  hash_itr itr;
  hash_reset(dst);
  hash_makeroom(dst, hash_card(src));
  for (hash_first(src, itr); hash_good(itr); hash_next(itr))
    {
      void *key = hash_key(itr);
      int val = hash_value(itr);
      hashkey_t k = itr.s->elts[itr.i].hkey;
      int *valp = _hash_mkfind_k(dst, key, k);
      *valp = val;
    }
}


size_t hash_find(hashtab *s, void *e, hashkey_t k)
{
  size_t index = k % hash_tabsz(s);
  size_t i = s->table[index];

  cmp_t cmp = hash_cmp(s);
  struct _hashelt *elts = s->elts;

  while (i != _S_END && cmp(e, elts[i].key))
    i = elts[i].next;

  return i;
}

void hash_makeroom(hashtab *s, size_t sz)
{
  if (USE_FAC * sz > hash_tabsz(s))
    {
      size_t *newtab;
      size_t newsz;
      size_t index, i, next;
      struct _hashelt *elts = s->elts;
      size_t *oldtab = s->table;

      newsz = 2 * USE_FAC * sz + 1;
      if (newsz % 3 == 0)
	newsz += 2;
      if (newsz % 5 == 0)
	newsz += 6;
      if (newsz % 7 == 0)
	newsz += 30;
      newtab = amalloc(newsz * sizeof(size_t));
      
      for (index = 0; index < newsz; index++)
	newtab[index] = _S_END;
      
      for (index = 0; index < hash_tabsz(s); index++)
	for (i = oldtab[index]; i != _S_END; i = next)
	  {
	    size_t ind = elts[i].hkey % newsz;
	    next = elts[i].next;
	    elts[i].next = newtab[ind];
	    newtab[ind] = i;
	  }
      
      s->table = newtab;
      hash_tabsz(s) = newsz;
      afree(oldtab);
    }
  
  if (sz > hash_eltsz(s))
    {
      size_t newsz = 2 * sz;
      size_t i;
      struct _hashelt *elts;
      elts = s->elts = arealloc(s->elts, newsz * sizeof(struct _hashelt));
      
      for (i = hash_eltsz(s); i < newsz; i++)
	elts[i].next = i + 1;
      
      elts[newsz - 1].next = s->free_elts;
      s->free_elts = hash_eltsz(s);
      
      hash_eltsz(s) = newsz;
    }
}


int hash_lookup(hashtab *s, void *key)
{
  hashkey_t k = hash_hash(s)(key);
  size_t i = hash_find(s, key, k);

  return (i == _S_END) ? 0 : s->elts[i].value;
}


int hash_insert(hashtab *ht, void *key, int value)
{
  int *valuep = hash_mkfind(ht, key);
  int oldvalue = *valuep;
  *valuep = value;
  return oldvalue;
}


int *_hash_mkfind_k(hashtab *ht, void *key, hashkey_t k)
{
  size_t i = hash_find(ht, key, k);
  
  if (i == _S_END)
    {
      size_t index;
      struct _hashelt *elts;
      
      hash_makeroom(ht, hash_card(ht) + 1);
      
      elts = ht->elts;
      index = k % hash_tabsz(ht);
      i = ht->free_elts;
      ht->free_elts = elts[i].next;
      
      elts[i].hkey = k;
      elts[i].key = key;
      elts[i].value = 0;
      
      elts[i].next = ht->table[index];
      ht->table[index] = i;
      
      hash_card(ht)++;
      
      hash_key_used = 1;
    }
  else
    {
      hash_key_used = (ht->elts[i].key == key);
    }
  return &(ht->elts[i].value);
}


int _hash_remove_k(hashtab *s, void *e, hashkey_t k)
{
  size_t index = k % hash_tabsz(s);
  size_t i = s->table[index];
  size_t prev = _S_END;
  
  cmp_t cmp = hash_cmp(s);
  struct _hashelt *elts = s->elts;
  
  while (i != _S_END && cmp(e, elts[i].key))
    {
      prev = i;
      i = elts[i].next;
    }
  
  if (i == _S_END)
    {
      hash_removed_key = NULL;
      return 0;
    }
  
  if (prev == _S_END)
    s->table[index] = elts[i].next;
  else
    elts[prev].next = elts[i].next;
  
  elts[i].next = s->free_elts;
  s->free_elts = i;
  
  hash_card(s)--;
  
  hash_removed_key = elts[i].key;

  return elts[i].value;
}


#if 0
list *hash_elemlist(hashtab *s)
{
  list *lst = l_newsz(hash_card(s));
  struct _hashelt *elts = hash_elts(s);
  size_t index;
  size_t i;

  for (index = 0; index < hash_tabsz(s); index++)
    for (i = s->table[index]; i != _S_END; i = elts[i].next)
      l_append(lst, elts[i].data);

  return lst;
}
#endif


void _hash_first(hashtab *s, hash_itr *itr)
{
  size_t index, i;

  itr->s = s;
  
  index = 0;
  while (index < hash_tabsz(s) && s->table[index] == _S_END)
    index++;

  if (index == hash_tabsz(s))
    {
      itr->i = _S_END;
      return;
    }
  
  i = s->table[index];
  itr->index = index;
  itr->i = i;
}

void _hash_next(hash_itr *itr)
{
  hashtab *s = itr->s;
  size_t index = itr->index;
  size_t i;
  
  index++;
  while (index < hash_tabsz(s) && s->table[index] == _S_END)
    index++;
  
  if (index == hash_tabsz(s))
    {
      itr->i = _S_END;
      return;
    }
  
  i = s->table[index];
  itr->index = index;
  itr->i = i;
}


void hash_print_stat(hashtab *s, size_t range)
{
  vector *stat = v_new_zero(range + 1);
  size_t index, i;
  
  for (index = 0; index < hash_tabsz(s); index++)
    {
      int count = 0;
      for (i = s->table[index]; i != _S_END; i = s->elts[i].next)
	count++;
      if (count > range)
	count = range;
      v_elem(stat, count)++;
    }
  
  puts("hash table distribution:");
  v_printnl(stat);
  
  v_free(stat);
}


void lincomb_add_multiple(hashtab *dst, int c, hashtab *lc, 
			  freekey_t freekey, copykey_t copykey)
{
  hash_itr itr;
  
  for (hash_first(lc, itr); hash_good(itr); hash_next(itr))
    {
      void *key = hash_key(itr);
      int value = hash_intvalue(itr);
      hashkey_t k = itr.s->elts[itr.i].hkey;
      
      int *valp = _hash_mkfind_k(dst, key, k);
      int newcoef = (*((int *) valp) += c * value);
      int hku = hash_key_used;
      
      if (newcoef == 0)
	{
	  _hash_remove_k(dst, key, k);
	  if (! hku)
	    freekey(hash_removed_key);
	  if (copykey == NULL)
	    freekey(key);
	}
      else
	{
	  if (copykey == NULL)
	    {
	      if (! hash_key_used)
		freekey(key);
	    }
	  else
	    {
	      if (hash_key_used)
		hash_key(itr) = copykey(key);
	    }
	}
    }
}

