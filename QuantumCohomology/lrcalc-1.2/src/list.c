/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include "list.h"
#include "claim.h"

list *l_newsz(size_t sz)
{
  list *lst = amalloc(sizeof(list));
  lst->array = amalloc(sizeof(void *) * sz);
  l_allocated(lst) = sz;
  l_length(lst) = 0;
  return lst;
}

void l_free(list *lst)
{
  afree(lst->array);
  afree(lst);
}

void l_makeroom(list *lst, size_t sz)
{
  if (l_allocated(lst) < sz)
    {
      l_allocated(lst) = 2 * sz;
      lst->array = arealloc(lst->array, l_allocated(lst) * sizeof(void *));
    }
}

void l_append(list *lst, void *e)
{
  l_makeroom(lst, l_length(lst) + 1);
  l_elem(lst, l_length(lst)++) = e;
}

void l_copy(list *dst, list *src)
{
  size_t i, len = l_length(src);
  l_makeroom(dst, len);
  l_length(dst) = len;
  for (i = 0; i < len; i++)
    l_elem(dst, i) = l_elem(src, i);
}

list *l_new_copy(list *lst)
{
  list *res = l_newsz(l_length(lst));
  l_copy(res, lst);
  return res;
}

void * l_deletelast(list *lst)
{
  void *x;
  size_t len = l_length(lst);
  claim(len > 0);
  x = l_elem(lst, len - 1);
  l_length(lst)--;
  return x;
}

void * l_delete(list *lst, size_t i)
{
  void *x = l_elem(lst, i);
  size_t j, len = l_length(lst);
  for (j = i; j < len - 1; j++)
    l_elem(lst, j) = l_elem(lst, j+1);
  l_length(lst)--;
  return x;
}

void * l_fastdelete(list *lst, size_t i)
{
  void *x = l_elem(lst, i);
  l_elem(lst, i) = l_elem(lst, l_length(lst) - 1);
  l_length(lst)--;
  return x;
}

void l_reset(list *lst)
{
  l_length(lst) = 0;
}
