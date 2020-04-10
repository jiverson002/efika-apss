/* SPDX-License-Identifier: MIT */
#include "efika/core.h"
#include "efika/impl.h"

#include "efika/core/blas.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"
#include "efika/impl/export.h"
#include "efika/impl/rename.h"


static int
bruteforce(val_t const minsim, Matrix const * const M, Matrix * const S)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* unpack /M/ */
  ind_t const m_nr = M->nr;
  ind_t const * const m_ia = M->ia;
  ind_t const * const m_ja = M->ja;
  val_t const * const m_a  = M->a;

  if (0 == m_nr)
    return 0;

  val_t * const spa = GC_calloc(m_nr, sizeof(*spa));

  //s_ia[0] = 0;
  for (ind_t i = 0/*, nnz = 0*/; i < m_nr - 1; i++) {
    BLAS_vsctr(m_ia[i + 1] - m_ia[i], m_a + m_ia[i], m_ja + m_ia[i], spa);

    for (ind_t j = i + 1; j < m_nr; j++) {
      val_t const d = BLAS_vdoti(m_ia[j + 1] - m_ia[j], m_a + m_ia[j],
                                 m_ja + m_ia[j], spa);
      if (d <= minsim) {
        //s_ja[nnz]  = j;
        //s_a[nnz++] = d;
      }
    }

    BLAS_vsctrz(m_ia[i + 1] - m_ia[i], m_ja + m_ia[i], spa);

    //s_ia[i + 1] = nnz;
  }

  //s_nnz = s_ia[s_nr];

  GC_free(spa);

  return 0;

  (void)S;
}


EFIKA_IMPL_EXPORT int
Impl_bruteforce(val_t const minsim, Matrix const * const M, Matrix * const S)
{
  if (!pp_all(M, S))
    return -1;

  return bruteforce(minsim, M, S);
}
