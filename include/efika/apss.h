/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_H
#define EFIKA_APSS_H 1

#include "efika/core.h"

#ifdef __cplusplus
extern "C" {
#endif

int EFIKA_apss_allpairs   (EFIKA_val_t const minsim,
                           EFIKA_Matrix * const M,
                           EFIKA_Matrix * const S);
int EFIKA_apss_allpairs_pp(EFIKA_val_t const minsim,
                           EFIKA_Matrix * const M);

int EFIKA_apss_bruteforce(EFIKA_val_t const minsim,
                          EFIKA_Matrix * const M,
                          EFIKA_Matrix * const S);

int EFIKA_apss_idxjoin   (EFIKA_val_t const minsim,
                          EFIKA_Matrix * const M,
                          EFIKA_Matrix * const S);
int EFIKA_apss_idxjoin_pp(EFIKA_val_t const minsim,
                          EFIKA_Matrix * const M);

#if 0
int EFIKA_apss_sfr(EFIKA_val_t const minsim,
                   EFIKA_Matrix const * const M,
                   EFIKA_Matrix const * const I,
                   EFIKA_Matrix       * const S);
#endif

#ifdef __cplusplus
}
#endif

#endif /* EFIKA_APSS_H */
