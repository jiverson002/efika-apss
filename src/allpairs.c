/* SPDX-License-Identifier: MIT */
#include <stdbool.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/allpairs.h"
#include "efika/apss/export.h"
#include "efika/apss/khash.h"
#include "efika/apss/rename.h"
#include "efika/core/blas.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"

/*----------------------------------------------------------------------------*/
/* ... */
/*----------------------------------------------------------------------------*/
#define UNKNOWN ((ind_t)-1)

#define WITH_KHASH
KHASH_MAP_INIT_INT64(m64, ind_t)

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
/*----------------------------------------------------------------------------*/ static inline ind_t
generate(
  ind_t const                  i,
  ind_t const * const restrict ia,
  ind_t const * const restrict ja,
  val_t const * const restrict a,
  ind_t const * const restrict ia1,
  ind_t const * const restrict ja1,
  val_t const * const restrict a1,
  ind_t       * const restrict marker,
  khash_t(m64) * const restrict hmarker,
  struct cand * const restrict tmpcnd
)
{
  ind_t cnt = 0;

  /* iterate through each non-zero column, j, for this row */
  for (ind_t jj = ia[i + 1]; jj > ia[i]; jj--) {
    ind_t const j = ja[jj - 1];
    val_t const v = a[jj - 1];

    /* iterate through the rows that have indexed column j */
    for (ind_t kk = ia1[j]; kk < ia1[j + 1]; kk++) {
      ind_t const k = ja1[kk];
      val_t const w = a1[kk];

      /* ignore candidates that come after i */
      if (k >= i)
        break;

#if defined(WITH_KHASH)
      int absent;
      khint_t key = kh_put(m64, hmarker, k, &absent);

      if (absent) {
        /* initialize solution matrix entry */
        tmpcnd[cnt].ind = k;
        tmpcnd[cnt].sim = v * w;

        kh_value(hmarker, key) = cnt++;
      } else {
        tmpcnd[kh_value(hmarker, key)].sim += v * w;
      }
#else
      ind_t const m = marker[k];
      switch (m) {
        case UNKNOWN:
        /* initialize solution matrix entry */
        tmpcnd[cnt].ind = k;
        tmpcnd[cnt].sim = v * w;

        /* populate marker */
        marker[k] = cnt++;
        break;

        default:
        /* update partial dot product for candidate row */
        tmpcnd[m].sim += v * w;
      }
#endif
    }
  }

  return cnt;

  (void)marker;
  (void)hmarker;
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
  khash_t(m64) * const restrict hmarker,
  val_t       * const restrict tmpspa,
  struct cand * const restrict tmpcnd
)
{
  ind_t ncnt = 0;

  BLAS_vsctr(ia[i + 1] - ia[i], a + ia[i], ja + ia[i], tmpspa);

  for (ind_t j = 0; j < cnt; j++) {
    ind_t const k = tmpcnd[j].ind;
    val_t       s = tmpcnd[j].sim;

#if defined(WITH_KHASH)
    /* reset hash map */
    for (khint_t key = kh_begin(hmarker); key != kh_end(hmarker); key++)
      if (kh_exist(hmarker, key))
        kh_del(m64, hmarker, key);
#else
    /* reset markers to unknown value */
    marker[k] = UNKNOWN;
#endif

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

  (void)marker;
  (void)hmarker;
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
  val_t       * const tmpspa = GC_calloc(m_nc, sizeof(*tmpspa));
  struct cand * const tmpcnd = GC_malloc(m_nr * sizeof(*tmpcnd));
  khash_t(m64) * const hmarker = kh_init(m64);

  /* initialize marker with unchecked value */
  for (ind_t i = 0; i < m_nr; i++)
    marker[i] = UNKNOWN;

  /* find similar neighbors for each query vector */
  s_ia[0] = 0;
  for (ind_t i = 0; i < m_nr; i++) {
    /* generate candidate vectors */
    ind_t const cnt = generate(i, m_ia, m_ja, m_a, i_ia, i_ja, i_a, marker,
                               hmarker, tmpcnd);

    /* verify candidate vectors */
    ind_t const ncnt = verify(minsim, cnt, i, m_ia, m_ja, m_ka, m_a, marker,
                              hmarker, tmpspa, tmpcnd);

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
  GC_free(tmpspa);
  GC_free(tmpcnd);
  kh_destroy(m64, hmarker);

  return 0;
}