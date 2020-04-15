/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_H
#define EFIKA_APSS_H 1

#include "efika/core.h"

/*----------------------------------------------------------------------------*/
/*! Various counters. */
/*----------------------------------------------------------------------------*/
extern unsigned long EFIKA_apss_ncand;
extern unsigned long EFIKA_apss_nprun;
extern unsigned long EFIKA_apss_nvdot;
extern unsigned long EFIKA_apss_nsims;
extern unsigned long EFIKA_apss_nmacs1;
extern unsigned long EFIKA_apss_nmacs2;

#ifdef __cplusplus
extern "C" {
#endif

int EFIKA_apss_allpairs   (EFIKA_val_t const minsim,
                           EFIKA_Matrix * const M,
                           EFIKA_Matrix * const S);
int EFIKA_apss_allpairs_pp(EFIKA_val_t const minsim,
                           EFIKA_Matrix * const M);

int EFIKA_apss_bruteforce   (EFIKA_val_t const minsim,
                             EFIKA_Matrix * const M,
                             EFIKA_Matrix * const S);
int EFIKA_apss_bruteforce_pp(EFIKA_val_t const minsim,
                             EFIKA_Matrix * const M);

int EFIKA_apss_idxjoin   (EFIKA_val_t const minsim,
                          EFIKA_Matrix * const M,
                          EFIKA_Matrix * const S);
int EFIKA_apss_idxjoin_pp(EFIKA_val_t const minsim,
                          EFIKA_Matrix * const M);

int EFIKA_apss_l2ap   (EFIKA_val_t const minsim,
                       EFIKA_Matrix * const M,
                       EFIKA_Matrix * const S);
int EFIKA_apss_l2ap_pp(EFIKA_val_t const minsim,
                       EFIKA_Matrix * const M);

int EFIKA_apss_mmjoin   (EFIKA_val_t const minsim,
                         EFIKA_Matrix * const M,
                         EFIKA_Matrix * const S);
int EFIKA_apss_mmjoin_pp(EFIKA_val_t const minsim,
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
