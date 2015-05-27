#ifndef _SCHUBLIB_H
#define _SCHUBLIB_H

#include <hashtab.h>

hashtab *trans(vector *w, int vars, hashtab *res);
hashtab *monk(int i, hashtab *slc, int rank);
hashtab *mult_schubert(vector *ww1, vector *ww2, int rank);
hashtab *mult_poly_schubert(hashtab *poly, vector *perm, int rank);

int num_inversions(vector *w);
int perm_group(vector *w);

list *all_strings(vector *dimvec);
list *all_perms(int n);

vector *string2perm(vector *str);
vector *str2dimvec(vector *str);
vector *perm2string(vector *perm, vector *dimvec);

hashtab *mult_str_schubert(vector *str1, vector *str2);

#endif
