/* SPDX-License-Identifier: MIT */
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/export.h"
#include "efika/apss/rename.h"
#include "efika/core/blas.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"

EFIKA_APSS_EXPORT int
apss_idxjoin(val_t const minsim, Matrix * const M, Matrix * const S)
{
  /* ...garbage collected function... */
  GC_func_init();

  if (!pp_all(M, S))
    return -1;

  /* unpack /M/ */
  ind_t const m_nr = M->nr;
  ind_t const m_nc = M->nc;
  ind_t const * const m_ia = M->ia;
  ind_t const * const m_ja = M->ja;
  val_t const * const m_a  = M->a;

  if (0 == m_nr)
    return 0;

  /* allocate memory for solution */
  ind_t nnz = 0;
  ind_t cap = m_nr;
  ind_t * const s_ia = GC_malloc((m_nr + 1) * sizeof(*s_ia));
  ind_t * s_ja = GC_malloc(cap * sizeof(*s_ja));
  val_t * s_a  = GC_malloc(cap * sizeof(*s_a));

  /* allocate temporary memory */
  val_t * const spa = GC_calloc(m_nc, sizeof(*spa));

  s_ia[0] = 0;
  for (ind_t i = 0; i < m_nr - 1; s_ia[i++] = nnz) {
    BLAS_vsctr(m_ia[i + 1] - m_ia[i], m_a + m_ia[i], m_ja + m_ia[i], spa);

    for (ind_t j = i + 1; j < m_nr; j++) {
      val_t const s = BLAS_vdoti(m_ia[j + 1] - m_ia[j], m_a + m_ia[j],
                                 m_ja + m_ia[j], spa);
      if (s >= minsim) {
        /* _vector_ resize */
        if (nnz >= cap) {
          cap *= 2;
          s_ja = GC_realloc(s_ja, cap * sizeof(*s_ja));
          s_a  = GC_realloc(s_a, cap * sizeof(*s_a));
        }

        s_ja[nnz]  = j;
        s_a[nnz++] = s;
      }
    }

    BLAS_vsctrz(m_ia[i + 1] - m_ia[i], m_ja + m_ia[i], spa);
  }

  /* record values in /S/ */
  /*S->fmt  = 0;*/
  S->sort   = ASC;
  /*S->symm = 0;*/
  S->nr     = m_nr;
  S->nc     = m_nr;
  S->nnz    = nnz;
  /*S->ncon = 0;*/
  S->ia     = s_ia;
  S->ja     = GC_realloc(s_ja, nnz * sizeof(*s_ja));
  S->a      = GC_realloc(s_a, nnz * sizeof(s_a));

  GC_free(spa);

  return 0;
}
