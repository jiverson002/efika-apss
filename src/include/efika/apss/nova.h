/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_NOVA_H
#define EFIKA_APSS_NOVA_H 1

#include "efika/apss/rename.h"

/*----------------------------------------------------------------------------*/
/* ... */
/*----------------------------------------------------------------------------*/
#define min(a, b)     (a < b ? a : b)
#define min3(a, b, c) min(min(a, b), c)
#define max(a, b)     (a > b ? a : b)

/*----------------------------------------------------------------------------*/
/*! Preprocessed data storage. */
/*----------------------------------------------------------------------------*/
struct pp_payld
{
  Matrix I;
  ind_t *ka;
  val_t *l;
  val_t *rowmax;
  val_t *colmax;
  val_t *pfxmax;
};

#endif /* EFIKA_APSS_NOVA_H */
