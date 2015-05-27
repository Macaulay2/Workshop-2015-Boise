#ifndef _ALLOC_H
#define _ALLOC_H

#include <stdlib.h>
#include <setjmp.h>

/*  Programs using the lrcalc library should set lrcalc_panic_frame
 *  with setjmp(lrcalc_panic_frame).  The lrcalc library will call
 *  longjmp(lrcalc_panic_frame, 1) if an "out of memory" event occurs.
 */
extern jmp_buf lrcalc_panic_frame;

void *amalloc(size_t size);
void *acalloc(size_t num, size_t size);
void *arealloc(void *p, size_t size);

#ifdef DEBUG
#define DEBUG_MEMORY
#endif

#if defined(DEBUG_MEMORY) || defined(SAGE)
void afree(void *);
#else
#define afree(p)	(free(p))
#endif

#ifdef DEBUG_MEMORY
void mem_report();
#define memory_report	mem_report()
#else
#define memory_report
#endif

#endif
