/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include "set.h"
#include "vector.h"


set *s_new(cmp_t cmp, hash_t hsh)
{
  set *s = amalloc(sizeof(set));
  s_cmp(s) = cmp;
  s_hash(s) = hsh;

  s_tabsz(s) = INIT_TABSZ;
  s_table(s) = amalloc(INIT_TABSZ * sizeof(size_t));

  s->elts_sz = INIT_ELTSZ;
  s_elts(s) = amalloc(INIT_ELTSZ * sizeof(struct _element));

  s_reset(s);

  return s;
}

void s_free(set *s)
{
  afree(s_table(s));
  afree(s_elts(s));
  afree(s);
}

void s_reset(set *s)
{
  size_t i;
  s_card(s) = 0;
  for (i = 0; i < s_tabsz(s); i++)
    s_table(s)[i] = _S_END;
  for (i = 0; i < s_eltsz(s); i++)
    s->elts[i].next = i + 1;
  s->elts[s_eltsz(s) - 1].next = _S_END;
  s->free_elts = 0;
}


size_t s_find(set *s, void *e, hashkey_t key)
{
  size_t index = key % s_tabsz(s);
  size_t i = s->table[index];

  cmp_t cmp = s_cmp(s);
  struct _element *elts = s->elts;

  while (i != _S_END && cmp(e, elts[i].data))
    i = elts[i].next;

  return i;
}

void s_makeroom(set *s, size_t sz)
{
  if (USE_FAC * sz > s_tabsz(s))
    {
      size_t *newtab;
      size_t newsz;
      size_t index, i, next;
      struct _element *elts = s->elts;
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
      
      for (index = 0; index < s_tabsz(s); index++)
	for (i = oldtab[index]; i != _S_END; i = next)
	  {
	    size_t ind = elts[i].key % newsz;
	    next = elts[i].next;
	    elts[i].next = newtab[ind];
	    newtab[ind] = i;
	  }
      
      s->table = newtab;
      s_tabsz(s) = newsz;
      afree(oldtab);
    }
  
  if (sz > s_eltsz(s))
    {
      size_t newsz = 2 * sz;
      size_t i;
      struct _element *elts;
      elts = s->elts = arealloc(s->elts, newsz * sizeof(struct _element));

      for (i = s_eltsz(s); i < newsz; i++)
	elts[i].next = i + 1;
      
      elts[newsz - 1].next = s->free_elts;
      s->free_elts = s_eltsz(s);
      
      s_eltsz(s) = newsz;
    }
}

void * s_insert(set *s, void *e)
{
  hashkey_t key = s_hash(s)(e);
  size_t i = s_find(s, e, key);

  if (i == _S_END)
    {
      size_t index;
      struct _element *elts;
      
      s_makeroom(s, s_card(s) + 1);
      
      elts = s->elts;
      index = key % s_tabsz(s);
      i = s->free_elts;
      s->free_elts = elts[i].next;
      
      elts[i].data = e;
      elts[i].key = key;

      elts[i].next = s->table[index];
      s->table[index] = i;

      s_card(s)++;

      return e;
    }
  else
    {
      return s->elts[i].data;
    }
}


void * s_has(set *s, void *e)
{
  hashkey_t key = s_hash(s)(e);
  size_t i = s_find(s, e, key);

  return (i == _S_END) ? NULL : s->elts[i].data;
}


void * s_remove(set *s, void *e)
{
  hashkey_t key = s_hash(s)(e);
  size_t index = key % s_tabsz(s);
  size_t i = s->table[index];
  size_t prev = _S_END;

  cmp_t cmp = s_cmp(s);
  struct _element *elts = s->elts;

  while (i != _S_END && cmp(e, elts[i].data))
    {
      prev = i;
      i = elts[i].next;
    }

  if (i == _S_END)
    return NULL;

  if (prev == _S_END)
    s->table[index] = elts[i].next;
  else
    elts[prev].next = elts[i].next;

  elts[i].next = s->free_elts;
  s->free_elts = i;

  s_card(s)--;

  return elts[i].data;
}


list *s_elemlist(set *s)
{
  list *lst = l_newsz(s_card(s));
  struct _element *elts = s_elts(s);
  size_t index;
  size_t i;

  for (index = 0; index < s_tabsz(s); index++)
    for (i = s->table[index]; i != _S_END; i = elts[i].next)
      l_append(lst, elts[i].data);

  return lst;
}


void *_s_first(set *s, set_itr *itr)
{
  size_t index, i;

  itr->s = s;
  
  index = 0;
  while (index < s_tabsz(s) && s->table[index] == _S_END)
    index++;

  if (index == s_tabsz(s))
    return NULL;

  i = s->table[index];
  itr->index = index;
  itr->i = i;
  return s->elts[i].data;
}

void *_s_next(set_itr *itr)
{
  set *s = itr->s;
  size_t index = itr->index;
  size_t i;
  
  index++;
  while (index < s_tabsz(s) && s->table[index] == _S_END)
    index++;
  
  if (index == s_tabsz(s))
    return NULL;

  i = s->table[index];
  itr->index = index;
  itr->i = i;
  return s->elts[i].data;
}


void s_print_stat(set *s, size_t range)
{
  vector *stat = v_new_zero(range + 1);
  size_t index, i;
  
  for (index = 0; index < s_tabsz(s); index++)
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
