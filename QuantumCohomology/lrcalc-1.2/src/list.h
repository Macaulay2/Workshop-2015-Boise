#ifndef _LIST_H
#define _LIST_H

#include "alloc.h"

typedef struct _list {
  void **array;
  size_t allocated;
  size_t length;
} list;

#define START_LENGTH 10

#define l_elem(lst, i)		((lst)->array[i])
#define l_length(lst)		((lst)->length)
#define l_allocated(lst)        ((lst)->allocated)


list *l_newsz(size_t sz);
#define l_new()		(l_newsz(10))
void l_free(list *lst);

void l_reset(list *lst);
void l_makeroom(list *lst, size_t size);

void l_append(list *lst, void *x);
void *l_deletelast(list *lst);
void l_insert(list *lst, size_t i, void *x);
void *l_delete(list *lst, size_t i);
void *l_fastdelete(list *lst, size_t i);

#define l_intelem(lst,i) ((int) (lst)->array[i])
#define l_appendint(lst,a) l_append(lst,(void *)(a))

void l_appendlist(list * lst1, list * lst2);
void l_copy(list * lst1, list * lst2);
list *l_new_copy(list * lst);

#endif
