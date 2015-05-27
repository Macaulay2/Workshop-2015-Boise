#ifndef _HASHTAB_H
#define _HASHTAB_H

#include <stdio.h>
#include <stdlib.h>

#ifndef _HASHKEY_T
#define _HASHKEY_T
typedef unsigned long hashkey_t;  /* should be 32 bit */
typedef int (*cmp_t) (void *, void *);  /* return 0 if equal */
typedef hashkey_t (*hash_t) (void *);
#endif

typedef struct {
  int card;
  cmp_t cmp;
  hash_t hash;
  
  size_t table_sz;
  size_t *table;
  
  size_t elts_sz;
  struct _hashelt *elts;
  size_t free_elts;
} hashtab;

struct _hashelt {
  size_t next;
  hashkey_t hkey;
  void *key;
  int value;
};

#define _S_END	((size_t) -1)

#define USE_FAC		2
#define INIT_TABSZ	2003
#define INIT_ELTSZ	100


#define hash_card(s)	((s)->card)
#define hash_cmp(s)	((s)->cmp)
#define hash_hash(s)	((s)->hash)

#define hash_tabsz(s)	((s)->table_sz)
#define hash_table(s)	((s)->table)
#define hash_eltsz(s)	((s)->elts_sz)
#define hash_elts(s)	((s)->elts)


extern int hash_key_used;
extern void *hash_removed_key;


hashtab *hash_new(cmp_t cm, hash_t hsh);
void hash_free(hashtab *ht);
void hash_reset(hashtab *ht);

void hash_copy(hashtab *ht1, hashtab *ht2);	/* ht1 := ht2 */
hashtab *hash_new_copy(hashtab *ht);


/* For int type values */

/* This macro renames hash_insert to avoids a name clash leading to a
   segmentation fault on Open Solaris.
   See also http://trac.sagemath.org/sage_trac/ticket/11563 */
#define hash_insert lrcalc_hash_insert

#define hash_lookupint(ht, key) \
  (hash_lookup((ht), (key)))
#define hash_insertint(ht, key, value)			\
  (hash_insert(ht, (key), (value)))
#define hash_mkfindint(ht, key) \
  (hash_mkfind((ht), (key)))


/* Returns value associated with key, or NULL if key is not in table. */

int hash_lookup(hashtab *ht, void *key);


/* Associates value to key.  The old value of key is returned.
 * If key is already in the hash table, the old key pointer is kept, 
 * and the given key pointer is not stored.  If so the global variable 
 * hash_key_used is set to 0, otherwise it is set to 1.  However, if
 * the old key pointer is physically equal to the given, hash_key_used is 
 * set to 1.
 *
 * Example:
 * 
 * oldvalue = hash_insert(ht, key, value);
 * if (! hash_key_used)
 *   free(key);
 * if (oldvalue != NULL && oldvalue != value)
 *   free(oldvalue);
 */

int hash_insert(hashtab *ht, void *key, int value);


/* Creates an entry in the hashtable with the given key.
 * A pointer to the value pointer of the hash table is returned.
 * hash_key_used is set like in hash_insert().
 *
 * Warning: The returned pointer may become invalid if more insertions
 * are made.
 *
 * Example:
 *
 * valuep = hash_mkfind(ht, key)
 * if (*valuep == NULL)
 *   {
 *     <calculate value>
 *     *valuep = value;
 *     if (! hash_key_used)
 *       free(key);
 *   }
 * value = *valuep;
 * <use value> */

#define hash_mkfind(ht, key)	(_hash_mkfind_k((ht), (key), hash_hash(ht)(key)))


/* Removes any entry with the given key from the hash table.
 * If such an entry exists, the key and value pointers for that entry 
 * are stored in the global variables hash_removed_key and 
 * hash_removed_value.  Otherwise these variables are set to NULL.
 */

#define hash_remove(ht, key)	(_hash_remove_k((ht), (key), hash_hash(ht)(e)))



#include "list.h"
list *hash_elemlist(hashtab *ht);

void hash_print_stat(hashtab *ht, size_t range);

typedef void (*freekey_t)(void *);
typedef void *(*copykey_t)(void *);
void lincomb_add_multiple(hashtab *res, int c, hashtab *lc, 
			  freekey_t freekey, copykey_t copykey);


typedef struct {
  hashtab *s;
  size_t index;
  size_t i;
} hash_itr;

#define hash_good(itr)	((itr).i != _S_END)
#define hash_key(itr)	((itr).s->elts[(itr).i].key)
#define hash_value(itr) ((itr).s->elts[(itr).i].value)
#define hash_intvalue(itr) ((int)(long) (itr).s->elts[(itr).i].value)

#define hash_first(s,itr)	(_hash_first((s), &(itr)))

#define hash_next(itr) \
  ((((itr).i = (itr).s->elts[(itr).i].next) == _S_END) \
   ? _hash_next(&(itr)) : 0)

void _hash_first(hashtab *s, hash_itr *itr);
void _hash_next(hash_itr *itr);

int *_hash_mkfind_k(hashtab *ht, void *key, hashkey_t k);
int _hash_remove_k(hashtab *ht, void *key, hashkey_t k);


#endif
