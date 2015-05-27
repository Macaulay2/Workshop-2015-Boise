/*  Littlewood-Richardson Calculator
 *  Copyright (C) 1999- Anders S. Buch (asbuch at math rutgers edu)
 *  See the file LICENSE for license information.
 */

#include "hashtab.h"
#include "vector.h"
#include "lincomb.h"


int lincomb_add_element(hashtab *dst, int c, void *elt, freekey_t freekey)
{
  hashkey_t k;
  int *valp;
  int newcoef;
  int hku;
  
  if (c == 0)
    return 0;
  
  k = hash_hash(dst)(elt);
  valp = _hash_mkfind_k(dst, elt, k);
  newcoef = (*valp += c);
  hku = hash_key_used;
  
  if (newcoef == 0)
    {
      _hash_remove_k(dst, elt, k);
      freekey(hash_removed_key);
    }
  return hku;
}

