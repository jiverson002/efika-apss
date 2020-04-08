/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_IMPL_H
#define EFIKA_IMPL_H 1

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

#include "efika/core.h"

typedef struct solution {
  unsigned p0;
  unsigned p1;
  float d;
} Solution;

typedef struct vector {
  size_t size;
  size_t capacity;
  Solution *mem;
} Vector;

static inline Solution vector_get(Vector * const v, size_t i) {
  return v->mem[i];
}

static inline void vector_push_back(Vector * const v, Solution const s) {
  if (v->size == v->capacity) {
    void *vmem = realloc(v->mem, 2 * v->capacity * sizeof(*v->mem));
    assert(vmem);

    v->capacity *= 2;
    v->mem = (Solution*)vmem;
  }

  v->mem[v->size++] = s;
}

static inline Vector vector_new(void) {
  Solution *vmem = (Solution*)malloc(1024 * sizeof(Solution));
  assert(vmem);
#ifdef __cplusplus
  return { 0, 1024, vmem };
#else
  return (Vector){ 0, 1024, vmem };
#endif
}

static inline void vector_delete(Vector * const v) {
  free(v->mem);
}

#ifdef __cplusplus
extern "C" {
#endif

int EFIKA_Impl_sfr0d(EFIKA_val_t const minsim,
                     EFIKA_Matrix const * const M,
                     Vector * const A);
int EFIKA_Impl_sfrkd(EFIKA_val_t const minsim,
                     EFIKA_Matrix const * const M,
                     EFIKA_Matrix const * const I,
                     Vector * const A);

#ifdef __cplusplus
}
#endif

#endif /* EFIKA_IMPL_H */
