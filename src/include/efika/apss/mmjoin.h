/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_MMJOIN_H
#define EFIKA_APSS_MMJOIN_H 1

#include "efika/apss/rename.h"

/*----------------------------------------------------------------------------*/
/* ... */
/*----------------------------------------------------------------------------*/
#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)

/*----------------------------------------------------------------------------*/
/*! Preprocessed data storage. */
/*----------------------------------------------------------------------------*/
struct pp_payld
{
  Matrix I;
  ind_t *ka;
  val_t *colmax;
  val_t *pfxmax;
};

#endif /* EFIKA_APSS_MMJOIN_H */
