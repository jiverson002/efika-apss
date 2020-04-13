/* SPDX-License-Identifier: MIT */
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/idxjoin.h"
#include "efika/apss/export.h"
#include "efika/apss/rename.h"
#include "efika/core/gc.h"

/*----------------------------------------------------------------------------*/
/*! Function to free any memory allocated during preprocessing. */
/*----------------------------------------------------------------------------*/
static void
pp_free(
  void * const ptr
)
{
  struct pp_payld * pp = (struct pp_payld *)ptr;
  Matrix_free(&(pp->I));
}

/*----------------------------------------------------------------------------*/
/*! Function to preprocess matrix for APSS. */
/*----------------------------------------------------------------------------*/
EFIKA_APSS_EXPORT int
apss_idxjoin_pp(
  val_t const minsim,
  Matrix * const M
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* ... */
  struct pp_payld * pp = GC_malloc(sizeof(*pp));

  /* unpack /pp/ */
  Matrix * const I = &(pp->I);

  /* init /I/ */
  int err = Matrix_init(I);
  GC_assert(!err);

  /* ... */
  err = Matrix_iidx(M, I);
  GC_assert(!err);

  /* record payload in /M/ */
  M->pp = pp;
  M->pp_free = &pp_free;

  return 0;

  (void)minsim;
}
