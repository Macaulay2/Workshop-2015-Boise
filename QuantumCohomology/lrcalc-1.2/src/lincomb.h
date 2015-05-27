#ifndef _LINCOMB_H
#define _LINCOMB_H

#if 0
typedef void (*freekey_t)(void *);
typedef void *(*copykey_t)(void *);
#endif
void lincomb_add_multiple(hashtab *res, int c, hashtab *lc, 
			  freekey_t freekey, copykey_t copykey);

/* Returns non-zero if elt is used: */
int lincomb_add_element(hashtab *dst, int c, void *elt, freekey_t freekey);

/* Number of non-zero coefficients */
int lincomb_card(hashtab *lc);

/* Return non-zero if linear combinations are different */
int lincomb_compare(hashtab *lc1, hashtab *lc2);

#endif
