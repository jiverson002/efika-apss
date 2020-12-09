/* SPDX-License-Identifier: MIT */
#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/bruteforce.h"
#include "efika/apss/rename.h"
#include "efika/core/export.h"
#include "efika/core/pp.h"

/*----------------------------------------------------------------------------*/
/*! Function to preprocess matrix for APSS. */
/*----------------------------------------------------------------------------*/
EFIKA_EXPORT int
apss_bruteforce_pp(
  val_t const minsim,
  Matrix * const M
)
{
  if (!M)
    return -1;

  return 0;

  (void)minsim;
  (void)M;
}
