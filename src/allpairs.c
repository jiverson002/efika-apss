/* SPDX-License-Identifier: MIT */
#include <stdbool.h>
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/allpairs.h"
#include "efika/apss/export.h"
#include "efika/apss/rename.h"
#include "efika/core/blas.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"

/*----------------------------------------------------------------------------*/
/* ... */
/*----------------------------------------------------------------------------*/
#define UNKNOWN ((ind_t)-1)

/*----------------------------------------------------------------------------*/
/*! Candidate object. */
/*----------------------------------------------------------------------------*/
struct cand
{
  /* solution matrix entry for candidate row */
  ind_t ind;
  val_t sim;
};

/*----------------------------------------------------------------------------*/
/*! */
/*----------------------------------------------------------------------------*/
static inline ind_t
generate(
  val_t const                  minsim,
  ind_t const                  i,
  ind_t const * const restrict ia,
  ind_t const * const restrict ja,
  val_t const * const restrict a,
  ind_t const * const restrict ia1,
  ind_t const * const restrict ja1,
  ind_t       * const restrict ra1,
  val_t const * const restrict a1,
  ind_t       * const restrict marker,
  val_t       * const restrict tmpmax,
  struct cand * const restrict tmpcnd
)
{
  ind_t cnt = 0;
  val_t mx = 0.0, rs1 = 0.0, sz1 = 0.0;

  /* precompute the minimum size */
  for (ind_t jj = ia[i]; jj < ia[i + 1]; jj++)
    mx = max(mx, a[jj]);
  sz1 = minsim / mx;

  /* precompute maximum dot product */
  for (ind_t jj = ia[i]; jj < ia[i + 1]; jj++)
    rs1 += a[jj] * tmpmax[ja[jj]];

  /* iterate through each non-zero column, j, for this row */
  for (ind_t jj = ia[i + 1]; jj > ia[i]; jj--) {
    ind_t const j = ja[jj - 1];
    val_t const v = a[jj - 1];

    /* ... */
    bool const allow_unknown = rs1 >= minsim;

    /* remove rows from index that are too short */
    for (; ra1[j] < ia1[j + 1]; ra1[j]++)
      if (ia[ja1[ra1[j]] + 1] - ia[ja1[ra1[j]]] >= sz1)
        break;

    /* iterate through the rows, k,  that have indexed column j */
    for (ind_t kk = ra1[j]; kk < ia1[j + 1]; kk++) {
      ind_t const k = ja1[kk];
      val_t const w = a1[kk];
      ind_t const m = marker[k];

      /* ignore candidates that come after i */
      if (k >= i)
        break;

      switch (m) {
        case UNKNOWN:
        if (allow_unknown) {
          /* initialize solution matrix entry */
          tmpcnd[cnt].ind = k;
          tmpcnd[cnt].sim = v * w;

          /* populate marker */
          marker[k] = cnt++;
        }
        break;

        default:
        /* update partial dot product for candidate row */
        tmpcnd[m].sim += v * w;
      }
    }

    /* ... */
    rs1 -= v * tmpmax[j];
  }

  return cnt;
}

/*----------------------------------------------------------------------------*/
/*! */
/*----------------------------------------------------------------------------*/
static inline ind_t
verify(
  val_t const                  minsim,
  ind_t const                  cnt,
  ind_t const                  i,
  ind_t const * const restrict ia,
  ind_t const * const restrict ja,
  ind_t const * const restrict ka,
  val_t const * const restrict a,
  ind_t       * const restrict marker,
  val_t       * const restrict tmpspa,
  struct cand * const restrict tmpcnd
)
{
  ind_t ncnt = 0;

  BLAS_vsctr(ia[i + 1] - ia[i], a + ia[i], ja + ia[i], tmpspa);

  for (ind_t j = 0; j < cnt; j++) {
    ind_t const k = tmpcnd[j].ind;
    val_t       s = tmpcnd[j].sim;

    /* reset markers to unknown value */
    marker[k] = UNKNOWN;

    /* compute the rest of the dot-product */
    s += BLAS_vdoti(ka[k] - ia[k], a + ia[k], ja + ia[k], tmpspa);

    /* retain only the cands that are at least minsim */
    if (s >= minsim) {
      tmpcnd[ncnt].ind   = k;
      tmpcnd[ncnt++].sim = s;
    }
  }

  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmpspa);

  return ncnt;
}

/*----------------------------------------------------------------------------*/
/*! Function to compute APSS of a matrix. */
/*----------------------------------------------------------------------------*/
EFIKA_APSS_EXPORT int
apss_allpairs(
  val_t const minsim,
  Matrix * const M,
  Matrix * const S
)
{
  /* ...garbage collected function... */
  GC_func_init();

  bool freepp = false;

  /* get /pp/ */
  if (!(M->pp)) {
    freepp = true;
    int err = apss_allpairs_pp(minsim, M);
    GC_assert(!err);
  }

  /* unpack /M/ */
  ind_t const m_nr = M->nr;
  ind_t const m_nc = M->nc;
  ind_t const * const m_ia = M->ia;
  ind_t const * const m_ja = M->ja;
  val_t const * const m_a  = M->a;
  struct pp_payld const * const pp = M->pp;

  /* unpack /pp/ */
  Matrix const * const I    = &(pp->I);
  ind_t  const * const m_ka = pp->ka;

  /* unpack /I/ */
  ind_t const         i_nr = I->nr;
  ind_t const * const i_ia = I->ia;
  ind_t const * const i_ja = I->ja;
  val_t const * const i_a  = I->a;

  /* allocate memory for solution */
  ind_t nnz = 0;
  ind_t cap = m_nr;
  ind_t * const s_ia = GC_malloc((m_nr + 1) * sizeof(*s_ia));
  ind_t * s_ja = GC_malloc(cap * sizeof(*s_ja));
  val_t * s_a  = GC_malloc(cap * sizeof(*s_a));

  /* allocate scratch memory */
  ind_t       * const marker = GC_malloc(m_nr * sizeof(*marker));
  val_t       * const tmpmax = GC_calloc(m_nc, sizeof(*tmpmax));
  val_t       * const tmpspa = GC_calloc(m_nc, sizeof(*tmpspa));
  struct cand * const tmpcnd = GC_malloc(m_nr * sizeof(*tmpcnd));
  ind_t       * const i_ra   = GC_malloc(i_nr * sizeof(*i_ra));

  /* initialize marker with unchecked value */
  for (ind_t i = 0; i < m_nr; i++)
    marker[i] = UNKNOWN;

  /* compute column maximums */
  for (ind_t i = 0; i < m_nr; i++)
    for(ind_t j = m_ia[i]; j < m_ia[i + 1]; j++)
      if (m_a[j] > tmpmax[m_ja[j]])
        tmpmax[m_ja[j]] = m_a[j];

  /* initialize ra1 */
  memcpy(i_ra, i_ia, i_nr * sizeof(*i_ra));

  /* find similar neighbors for each query vector */
  s_ia[0] = 0;
  for (ind_t i = 0; i < m_nr; i++) {
    /* generate candidate vectors */
    ind_t const cnt = generate(minsim, i, m_ia, m_ja, m_a, i_ia, i_ja, i_ra,
                               i_a, marker, tmpmax, tmpcnd);

    /* verify candidate vectors */
    ind_t const ncnt = verify(minsim, cnt, i, m_ia, m_ja, m_ka, m_a, marker,
                              tmpspa, tmpcnd);

    /* update successor column max */
    for (ind_t j = m_ia[i]; j < m_ia[i + 1]; j++)
      tmpmax[m_ja[j]] = max(m_a[j], tmpmax[m_ja[j]]);

    /* _vector_ resize */
    if (nnz + ncnt >= cap) {
      do {
        cap *= 2;
      } while (nnz + ncnt >= cap);
      s_ja = GC_realloc(s_ja, cap * sizeof(*s_ja));
      s_a  = GC_realloc(s_a, cap * sizeof(*s_a));
    }
    /* populate solution index */
    for (ind_t j = 0; j < ncnt; j++) {
      s_ja[nnz]  = tmpcnd[j].ind;
      s_a[nnz++] = tmpcnd[j].sim;
    }
    s_ia[i + 1] = nnz;
  }

  /* record values in /S/ */
  /*S->fmt  = 0;*/
  S->sort   = NONE;
  /*S->symm = 0;*/
  S->nr     = m_nr;
  S->nc     = m_nr;
  S->nnz    = nnz;
  /*S->ncon = 0;*/
  S->ia     = s_ia;
  S->ja     = GC_realloc(s_ja, nnz * sizeof(*s_ja));
  S->a      = GC_realloc(s_a, nnz * sizeof(s_a));

  /* free preprocessed data (if it was allocated here) */
  if (freepp) {
    M->pp_free(M->pp);
    M->pp = NULL;
  }

  /* free scratch memory */
  GC_free(marker);
  GC_free(tmpmax);
  GC_free(tmpspa);
  GC_free(tmpcnd);
  GC_free(i_ra);

  return 0;
}
