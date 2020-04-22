/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_NOVA_H
#define EFIKA_APSS_NOVA_H 1

#include "efika/core.h"

#include "efika/apss/rename.h"

/*----------------------------------------------------------------------------*/
/* ... */
/*----------------------------------------------------------------------------*/
#define min(a, b)        ((a) < (b) ? (a) : (b))
#define min3(a, b, c)    min(min(a, b), c)
#define min4(a, b, c, d) min(min(a, b), min(c, d))
#define max(a, b)        ((a) > (b) ? (a) : (b))
#define max3(a, b, c)    max(max(a, b), c)
#define sum(a, b)        ((a) + (b))

/*----------------------------------------------------------------------------*/
/*! Preprocessed data storage. */
/*----------------------------------------------------------------------------*/
struct pp_payld
{
  Matrix I;
  ind_t *ka;
  ind_t *ra;

  val_t *m_max;
  val_t *m_sum;
  val_t *m_sqr;
  val_t *m_len;
  val_t *m_rs1;

  val_t *i_max;
  val_t *i_sum;
  val_t *i_sqr;
  val_t *i_len;
  val_t *i_rs1;

  ind_t *psplit;
  val_t *pscore;
};

#endif /* EFIKA_APSS_NOVA_H */
