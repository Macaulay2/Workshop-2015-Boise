#ifndef _CLAIM_H_
#define _CLAIM_H_

/* The macro claim works like the standard C macro assert except that when
   it fails it provokes a segmentation fault, which makes it possible to 
   use a gdb */

#ifdef DEBUG
#define CLAIM
#endif

#include <stdio.h>
#ifdef CLAIM
#define claim(ex) ((void) ((ex) || (fprintf(stderr, "%s:%d: Claim `%s' failed\n", __FILE__, __LINE__, #ex), (*(char *)-1)=7)))
#else
#define claim(ex) ((void) 0)
#endif

#endif
