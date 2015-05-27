/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

/*  This module replaces alloc.c when lrcalc is linked programs that
 *  require lrcalc functions to free up memory and return in case of
 *  an out-of-memory event.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "alloc.h"
jmp_buf lrcalc_panic_frame;

typedef struct mlink {
  struct mlink *next;
  struct mlink *prev;
} mlink;

mlink root;

void panic()
{
  mlink *link, *link0;
  link = root.next;
  while (link)
    {
      link0 = link;
      link = link->next;
      free(link0);
    }
  root.next = NULL;
  root.prev = NULL;
  longjmp(lrcalc_panic_frame, 1);
}

void *amalloc(size_t size)
{
  mlink *p = malloc(size + sizeof(mlink));
  if (p == NULL)
    panic();
  p->next = root.next;
  p->prev = &root;
  root.next = p;
  if (p->next)
    p->next->prev = p;
  return ((void *)p) + sizeof(mlink);
}

void *acalloc(size_t size, size_t num)
{
  mlink *p = calloc(size*num + sizeof(mlink), 1);
  if (p == NULL)
    panic();
  p->next = root.next;
  p->prev = &root;
  root.next = p;
  if (p->next)
    p->next->prev = p;
  return ((void *)p) + sizeof(mlink);
}

void *arealloc(void *pp, size_t size)
{
  mlink *p = (mlink *) (pp - sizeof(mlink));
  p->prev->next = p->next;
  if (p->next)
    p->next->prev = p->prev;
  p = realloc(p, size + sizeof(mlink));
  if (p == NULL)
    {
      free(p);
      panic();
    }
  p->next = root.next;
  p->prev = &root;
  root.next = p;
  if (p->next)
    p->next->prev = p;
  return ((void *)p) + sizeof(mlink);
}

void afree(void *pp)
{
  mlink *p = (mlink *) (pp - sizeof(mlink));
  p->prev->next = p->next;
  if (p->next)
    p->next->prev = p->prev;
  free(p);
}

