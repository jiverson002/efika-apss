/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_IMPL_H
#define EFIKA_IMPL_H 1

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

typedef struct problem {
  float t;
  unsigned k;
  unsigned n;
  float (*mem)[5];
} Problem;

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
  return (Vector){ 0, 1024, vmem };
}

static inline void vector_delete(Vector * const v) {
  free(v->mem);
}

#ifdef __cplusplus
extern "C" {
#endif

void EFIKA_Impl_sfr1d(Problem const P, Vector * const A);
void EFIKA_Impl_sfrkd(Problem const P, Vector * const A);

#ifdef __cplusplus
}
#endif

#endif /* EFIKA_IMPL_H */
