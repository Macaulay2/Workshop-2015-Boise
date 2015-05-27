#ifndef _SYMFCN_H
#define _SYMFCN_H

#include <hashtab.h>
#include <vector.h>

int part_itr_sz(vector *part);
int part_itr_sub(vector *part, vector *outer);
int part_itr_between(vector *part, vector *inner, vector *outer);

int part_length(vector *p);
vector *part_conjugate(vector *p);
int part_subset(vector *p1, vector *p2);

long long lrcoef(vector *outer, vector *inner1, vector *inner2);

hashtab *skew(vector *outer, vector *inner, int maxrows);
hashtab *mult(vector *sh1, vector *sh2, int maxrows);
hashtab *coprod(vector *part, int all);

hashtab *schur_lc_mult(hashtab *lc1, hashtab *lc2, int maxrows);

list *quantum_reduce(hashtab* s, int rows, int cols);
void fusion_reduce(hashtab *lc, int rows, int cols, int opt_zero);


typedef struct {
  vector *outer;
  vector *inner;
  vector *conts;
  int maxrows;
  vector *conjugate;
  int rows;
  int cols;
  int matrix[1];
} skewtab;

skewtab *st_new(vector *outer, vector *inner, vector *conts, int maxrows);
int st_next(skewtab *st);
void st_print(skewtab *st);
void st_free(skewtab *st);

#endif
