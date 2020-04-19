/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_BIT_H
#define EFIKA_APSS_BIT_H 1

#include "efika/core.h"

#include "efika/core/rename.h"

typedef val_t (*BIT_op)(val_t, val_t);

static inline ind_t
BIT_lsb(ind_t const i) {
  return i & -i;
}

static inline val_t
BIT_sum(ind_t i, val_t const * const restrict a, BIT_op const op) {
  val_t res = 0.0;
  for (i++; i > 0; i -= BIT_lsb(i))
    res = op(a[i], res);
  return res;
}

static inline void
BIT_add(val_t k, ind_t n, ind_t i, val_t * const a, BIT_op const op) {
  for (i++; i <= n; i += BIT_lsb(i))
    a[i] = op(k, a[i]);
}

static inline val_t
BIT_op_sum(val_t a, val_t b) {
  return a + b;
}

static inline val_t
BIT_op_max(val_t a, val_t b) {
  return a > b ? a : b;
}

#endif /* EFIKA_APSS_BIT_H */
