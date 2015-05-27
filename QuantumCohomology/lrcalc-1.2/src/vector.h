#ifndef _VECTOR_H
#define _VECTOR_H

#include <stdlib.h>
#include <string.h>

typedef struct {
  size_t length;
  int array[1];
} vector;

#define v_length(v)	((v)->length)
#define v_elem(v,i)	((v)->array[i])

vector *v_new(int length);
vector *v_new_init(int length, ...);
vector *v_new_zero(int length);
#define v_free(v)	(afree(v))
void v_free1(vector *v);  /* a compiled version */

vector *v_new_copy(vector *v);
#define v_copy(d, v)            \
        memcpy((d)->array, (v)->array, (d)->length * sizeof(int))

#define v_set_zero(v)	(memset((v)->array, 0, (v)->length * sizeof(int)))
int v_sum(vector *v);
int v_lesseq(vector *v1, vector *v2);
void v_mult(vector *dst, int c, vector *src);
void v_div(vector *dst, vector *src, int c);
int v_max(vector *v);

int v_elm_cmp(const void *, const void *);
#define v_sort(v)	\
  qsort((v)->array, (v)->length, sizeof(int), v_elm_cmp);
#define v_sort_part(v, first, n)	\
  qsort((v)->array + (first), (n), sizeof(int), v_elm_cmp);

int i_gcd(int x, int y);
int v_gcd(vector *v);

void v_print(vector *v);
void v_printnl(vector *v);

#include "set.h"
#include "hashtab.h"
int v_cmp(vector *v1, vector *v2);
hashkey_t v_hash(vector *v);

void print_vec_set(set *s);
void free_vec_set(set *s);

void print_vec_lincomb(hashtab *ht, int opt_zero);
void free_vec_lincomb(hashtab *ht);


typedef struct {
  vector *first, *second;
} vecpair;

#define vp_first(vp)	((vp)->first)
#define vp_second(vp)	((vp)->second)

vecpair *vp_new(vector *v1, vector *v2);
vecpair *vp_new_unordered(vector *v1, vector *v2);
void vp_free(vecpair *vp);

void print_vp_lincomb(hashtab *ht);
void free_vp_lincomb(hashtab *ht);

int vp_cmp(vecpair *vp1, vecpair *vp2);
hashkey_t vp_hash(vecpair *vp);


#include "list.h"
void print_vec_list(list *lst);
void free_vec_list(list *lst);

list *find_extreme_vectors(list *veclist, int take_max);

#endif
