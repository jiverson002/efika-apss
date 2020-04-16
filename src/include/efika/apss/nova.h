/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_NOVA_H
#define EFIKA_APSS_NOVA_H 1

#include "efika/apss/rename.h"

/*----------------------------------------------------------------------------*/
/* ... */
/*----------------------------------------------------------------------------*/
#define min(a, b)     ((a) < (b) ? (a) : (b))
#define min3(a, b, c) min(min(a, b), c)
#define max(a, b)     ((a) > (b) ? (a) : (b))

/*----------------------------------------------------------------------------*/
/*! Preprocessed data storage. */
/*----------------------------------------------------------------------------*/
struct pp_payld
{
  Matrix I;
  ind_t *ka;
  ind_t *ra;
  val_t *m_rs1;
  val_t *m_rs3;
  val_t *i_rs3;
  val_t *rowmax;
  val_t *pfxmax;
  val_t *pscore;
};

#endif /* EFIKA_APSS_NOVA_H */
