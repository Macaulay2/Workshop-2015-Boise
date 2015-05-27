#ifndef _SET_H
#define _SET_H

#include <stdio.h>
#include <stdlib.h>

#ifndef _HASHKEY_T
#define _HASHKEY_T
typedef unsigned long hashkey_t;  /* should be 32 bit */
typedef int (*cmp_t) (void *, void *);  /* return 0 if equal */
typedef hashkey_t (*hash_t) (void *);
#endif

typedef struct _set {
  int card;
  cmp_t cmp;
  hash_t hash;

  size_t table_sz;
  size_t *table;

  size_t elts_sz;
  struct _element *elts;
  size_t free_elts;
} set;

struct _element {
  size_t next;
  void *data;
  hashkey_t key;
};

#define _S_END	((size_t) -1)


#define USE_FAC		2
#define INIT_TABSZ	2003
#define INIT_ELTSZ	100


#define s_card(s)	((s)->card)
#define s_cmp(s)	((s)->cmp)
#define s_hash(s)	((s)->hash)

#define s_tabsz(s)	((s)->table_sz)
#define s_table(s)	((s)->table)
#define s_eltsz(s)	((s)->elts_sz)
#define s_elts(s)	((s)->elts)


set *s_new(cmp_t cm, hash_t hsh);
void s_free(set *s);
void s_reset(set *s);

void s_copy(set *s1, set *s2);	     /* s1 = s2 */
set *s_new_copy(set *s);

void *s_has(set *s, void *e);
void *s_insert(set *s, void *e);
void *s_remove(set *s, void *e);

void s_union(set *s1, set *s2);      /* s1 = s1 union s2 */
void s_intersect(set *s1, set *s2);  /* s1 = s1 intersect s2 */
void s_minus(set *s1, set *s2);      /* s1 = s1 - s2 */
void s_symdiff(set *s1, set *s2);    /* s1 = (s1 - s2) union (s2 - s1) */

#include "list.h"
list *s_elemlist(set *s);

void s_print_stat(set *s, size_t range);


typedef struct {
  set *s;
  size_t index;
  size_t i;
} set_itr;

#define s_first(s,itr)		(_s_first((s), &itr))

#define s_next(itr) \
  ((itr.s->elts[itr.i].next == _S_END) \
   ? _s_next(&itr) \
   : itr.s->elts[itr.i = itr.s->elts[itr.i].next].data)

void *_s_first(set *s, set_itr *itr);
void *_s_next(set_itr *itr);


#endif
