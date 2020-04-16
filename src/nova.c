/* SPDX-License-Identifier: MIT */
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "efika/apss.h"

#include "efika/apss/nova.h"
#include "efika/apss/export.h"
#include "efika/apss/rename.h"
#include "efika/core/blas.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"

/*----------------------------------------------------------------------------*/
/* ... */
/*----------------------------------------------------------------------------*/
#define UNKNOWN ((ind_t)-1)
#define PRUNED  ((ind_t)-2)

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
  ind_t const * const restrict m_ia,
  ind_t const * const restrict m_ja,
  ind_t const * const restrict m_ka,
  ind_t const * const restrict m_ra,
  val_t const * const restrict m_a,
  val_t const * const restrict m_rs1,
  val_t const * const restrict m_rs3,
  ind_t const * const restrict i_ia,
  ind_t const * const restrict i_ja,
  val_t const * const restrict i_a,
  val_t const * const restrict i_rs3,
  ind_t       * const restrict marker,
  val_t       * const restrict tmprs3,
  struct cand * const restrict tmpcnd
)
{
  ind_t cnt = 0;
  val_t rs1 = 1.0, rs3 = 1.0;

  /* iterate through each non-zero column, j, for this row */
  for (ind_t jj = m_ia[i + 1]; jj > m_ia[i]; jj--) {
    ind_t const j    = m_ja[jj - 1];
    val_t const v    = m_a[jj - 1];
    ind_t const i_ra = m_ra[jj - 1];

    /* determine if any new candidates should be allowed */
    bool const allow_unknown  = /* jj > m_ka[i] */ true
                             && min(rs1, rs3) >= minsim;

    /* retrieve precomputed values */
    rs1 = m_rs1[jj - 1];
    rs3 = m_rs3[jj - 1];

    /* ... */
    tmprs3[j] = rs3;

    /* iterate through the rows, k,  that have indexed column j */
    for (ind_t kk = i_ra; kk < i_ia[j + 1]; kk++) {
      ind_t const k = i_ja[kk];
      val_t const w = i_a[kk];
      ind_t const m = marker[k];

      /* ignore candidates that come after i */
      if (k >= i)
        break;

      switch (m) {
        case PRUNED:
        break;

        case UNKNOWN:
        if (allow_unknown) {
          /* initialize solution matrix entry */
          tmpcnd[cnt].ind = k;
          tmpcnd[cnt].sim = v * w;

          /* populate marker */
          marker[k] = cnt++;

          /* update global counter */
          apss_nmacs1++;
        }
        break;

        default:
        /* update partial dot product for candidate row */
        tmpcnd[m].sim += v * w;

        /* Anastasiu pruning */
        if (tmpcnd[m].sim + rs3 * i_rs3[kk] < minsim)
          marker[k] = PRUNED;

        /* update global counter */
        apss_nmacs1++;
      }
    }
  }

  return cnt;

  (void)m_ka;
  (void)i_rs3;
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
  val_t const * const restrict rs3,
  val_t const * const restrict rowmax,
  val_t const * const restrict pfxmax,
  val_t const * const restrict pscore,
  ind_t       * const restrict marker,
  val_t       * const restrict tmprs3,
  val_t       * const restrict tmpspa,
  struct cand * const restrict tmpcnd
)
{
  ind_t ncnt = 0, len = ia[i + 1] - ia[i];

  BLAS_vsctr(ia[i + 1] - ia[i], a + ia[i], ja + ia[i], tmpspa);

  for (ind_t j = 0; j < cnt; j++) {
    ind_t const k = tmpcnd[j].ind;
    val_t       s = tmpcnd[j].sim;
    ind_t const m = marker[k];

    /* reset markers to unknown value */
    marker[k] = UNKNOWN;

    /* ... */
    if (PRUNED == m) {
      /* update global counter */
      apss_nprun++;

      continue;
    }

    /* Bayardo filter */
    if (s + min(ka[k] - ia[k], len) * rowmax[i] * pfxmax[k] < minsim)
      continue;

    /* pscore filter */
    if (s + pscore[k] < minsim)
      continue;

    /* Anastasiu dot product */
    for(ind_t jj = ka[k]; jj > ia[k]; jj--){
      ind_t const jjj = jj - 1;
      if (tmpspa[ja[jjj]] > 0.0) {
        s += tmpspa[ja[jjj]] * a[jjj];

        /* update global counter */
        apss_nmacs2++;

        if (s + tmprs3[ja[jjj]] + rs3[jjj] < minsim) {
          s = 0.0;
          break;
        }
      }
    }

    /* retain only the cands that are at least minsim */
    if (s >= minsim) {
      tmpcnd[ncnt].ind   = k;
      tmpcnd[ncnt++].sim = s;
    }

    /* update global counter */
    apss_nvdot++;
  }

  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmprs3);
  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmpspa);

  return ncnt;
}

/*----------------------------------------------------------------------------*/
/*! Function to compute APSS of a matrix. */
/*----------------------------------------------------------------------------*/
EFIKA_APSS_EXPORT int
apss_nova(
  val_t const minsim,
  Matrix * const M,
  Matrix * const S
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* ... */
  if (!pp_all(M, S))
    return -1;

  bool freepp = false;

  /* get /pp/ */
  if (!(M->pp)) {
    freepp = true;
    int err = apss_nova_pp(minsim, M);
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
  Matrix const * const I = &(pp->I);
  ind_t  const * const m_ka = pp->ka;
  ind_t  const * const m_ra = pp->ra;
  val_t  const * const m_rs1 = pp->m_rs1;
  val_t  const * const m_rs3 = pp->m_rs3;
  val_t  const * const i_rs3 = pp->i_rs3;
  val_t  const * const rowmax = pp->rowmax;
  val_t  const * const pfxmax = pp->pfxmax;
  val_t  const * const pscore = pp->pscore;

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
  val_t       * const tmprs3 = GC_calloc(m_nc, sizeof(*tmprs3));
  val_t       * const tmpspa = GC_calloc(m_nc, sizeof(*tmpspa));
  struct cand * const tmpcnd = GC_malloc(m_nr * sizeof(*tmpcnd));

  /* initialize marker with unchecked value */
  for (ind_t i = 0; i < m_nr; i++)
    marker[i] = UNKNOWN;

  /* reset global counters */
  apss_ncand = 0;
  apss_nprun = 0;
  apss_nvdot = 0;
  apss_nsims = 0;
  apss_nmacs1 = 0;
  apss_nmacs2 = 0;

  /* find similar neighbors for each query vector */
  s_ia[0] = 0;
  for (ind_t i = 0; i < m_nr; i++) {
    /* generate candidate vectors */
    ind_t const cnt = generate(minsim, i, m_ia, m_ja, m_ka, m_ra, m_a, m_rs1,
                               m_rs3, i_ia, i_ja, i_a, i_rs3, marker, tmprs3,
                               tmpcnd);

    /* verify candidate vectors */
    ind_t const ncnt = verify(minsim, cnt, i, m_ia, m_ja, m_ka, m_a, m_rs3,
                              rowmax, pfxmax, pscore, marker, tmprs3, tmpspa,
                              tmpcnd);

    /* update global counters */
    apss_ncand += cnt;
    apss_nsims += ncnt;

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
  GC_free(tmprs3);
  GC_free(tmpspa);
  GC_free(tmpcnd);

  return 0;
}
