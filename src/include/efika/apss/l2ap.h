/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_L2AP_H
#define EFIKA_APSS_L2AP_H 1

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
  val_t *l;
  val_t *rowmax;
  val_t *pfxmax;
  val_t *pscore;
};

#endif /* EFIKA_APSS_L2AP_H */
