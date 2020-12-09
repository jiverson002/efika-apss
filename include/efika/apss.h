/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_H
#define EFIKA_APSS_H 1

#include "efika/core.h"

#include "efika/core/export.h"

/*----------------------------------------------------------------------------*/
/*! Various counters. */
/*----------------------------------------------------------------------------*/
EFIKA_EXPORT extern unsigned long EFIKA_apss_ncand;
EFIKA_EXPORT extern unsigned long EFIKA_apss_nprun;
EFIKA_EXPORT extern unsigned long EFIKA_apss_nvdot;
EFIKA_EXPORT extern unsigned long EFIKA_apss_nsims;
EFIKA_EXPORT extern unsigned long EFIKA_apss_nmacs1;
EFIKA_EXPORT extern unsigned long EFIKA_apss_nmacs2;

#ifdef __cplusplus
extern "C" {
#endif

EFIKA_EXPORT int EFIKA_apss_allpairs   (EFIKA_val_t const minsim,
                                        EFIKA_Matrix * const M,
                                        EFIKA_Matrix * const S);
EFIKA_EXPORT int EFIKA_apss_allpairs_pp(EFIKA_val_t const minsim,
                                        EFIKA_Matrix * const M);

EFIKA_EXPORT int EFIKA_apss_bruteforce   (EFIKA_val_t const minsim,
                                          EFIKA_Matrix * const M,
                                          EFIKA_Matrix * const S);
EFIKA_EXPORT int EFIKA_apss_bruteforce_pp(EFIKA_val_t const minsim,
                                          EFIKA_Matrix * const M);

EFIKA_EXPORT int EFIKA_apss_idxjoin   (EFIKA_val_t const minsim,
                                       EFIKA_Matrix * const M,
                                       EFIKA_Matrix * const S);
EFIKA_EXPORT int EFIKA_apss_idxjoin_pp(EFIKA_val_t const minsim,
                                       EFIKA_Matrix * const M);

EFIKA_EXPORT int EFIKA_apss_l2ap   (EFIKA_val_t const minsim,
                                    EFIKA_Matrix * const M,
                                    EFIKA_Matrix * const S);
EFIKA_EXPORT int EFIKA_apss_l2ap_pp(EFIKA_val_t const minsim,
                                    EFIKA_Matrix * const M);

EFIKA_EXPORT int EFIKA_apss_mmjoin   (EFIKA_val_t const minsim,
                                      EFIKA_Matrix * const M,
                                      EFIKA_Matrix * const S);
EFIKA_EXPORT int EFIKA_apss_mmjoin_pp(EFIKA_val_t const minsim,
                                      EFIKA_Matrix * const M);

EFIKA_EXPORT int EFIKA_apss_nova   (EFIKA_val_t const minsim,
                                    EFIKA_Matrix * const M,
                                    EFIKA_Matrix * const S);
EFIKA_EXPORT int EFIKA_apss_nova_pp(EFIKA_val_t const minsim,
                                    EFIKA_Matrix * const M);

#if 0
EFIKA_EXPORT int EFIKA_apss_sfr(EFIKA_val_t const minsim,
                                EFIKA_Matrix const * const M,
                                EFIKA_Matrix const * const I,
                                EFIKA_Matrix       * const S);
#endif

#ifdef __cplusplus
}
#endif

#endif /* EFIKA_APSS_H */
