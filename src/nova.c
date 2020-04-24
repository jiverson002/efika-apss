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
  ind_t const * const restrict m_ra,
  val_t const * const restrict m_a,
  val_t const * const restrict m_max,
  val_t const * const restrict m_sum,
  val_t const * const restrict m_len,
  ind_t const * const restrict i_ia,
  ind_t const * const restrict i_ja,
  val_t const * const restrict i_a,
  val_t const * const restrict i_max,
  val_t const * const restrict i_sum,
  val_t const * const restrict i_len,
  val_t const * const restrict i_rs1,
  ind_t const * const restrict psplit,
  ind_t       * const restrict marker,
  struct cand * const restrict tmpcnd
)
{
  ind_t cnt = 0;

  /* iterate through each non-zero column, j, for this row */
  for (ind_t jj = m_ia[i + 1]; jj > m_ia[i]; jj--) {
    ind_t const j = m_ja[jj - 1];
    val_t const v = m_a[jj - 1];

    /* retrieve precomputed values */
    ind_t const ra  = m_ra[jj - 1];
    val_t const max = m_max[jj - 1];
    val_t const sum = m_sum[jj - 1];
    val_t const len = m_len[jj - 1];

    /* determine if any new candidates should be allowed */
    bool const allow_unknown = jj > psplit[i];

    /* iterate through the rows, k,  that have indexed column j */
    for (ind_t kk = ra; kk < i_ia[j + 1]; kk++) {
      ind_t const k = i_ja[kk];
      val_t const w = i_a[kk];
      ind_t       m = marker[k];

      /* ignore candidates that come after i */
      if (k >= i)
        break;

      switch (m) {
        case PRUNED:
        break;

        case UNKNOWN:
        if (!allow_unknown)
          break;

        /* initialize solution matrix entry */
        tmpcnd[cnt].ind = k;
        tmpcnd[cnt].sim = 0.0;

        /* populate marker and updated m */
        m = marker[k] = cnt++;

        /* fall-through */

        default:
#if 1
        (void)max;
        (void)sum;
        (void)i_max;
        (void)i_sum;
        (void)i_rs1;
        if ( tmpcnd[m].sim + len * i_len[kk] < minsim ) /* Anastasiu */
#else
        if ( tmpcnd[m].sim + len * i_len[kk] < minsim   /* Anastasiu */
          || tmpcnd[m].sim + i_rs1[kk]       < minsim   /* Bayardo   */
          || tmpcnd[m].sim + max * i_sum[kk] < minsim   /* Awekar    */
          || tmpcnd[m].sim + i_max[kk] * sum < minsim ) /* ...       */
#endif
        {
          marker[k] = PRUNED;
          break;
        }

        /* update partial dot product for candidate row */
        tmpcnd[m].sim += v * w;

        /* update global counter */
        apss_nmacs1++;
      }
    }
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
  val_t const * const restrict max,
  val_t const * const restrict sum,
  val_t const * const restrict len,
  val_t const * const restrict rs1,
  val_t const * const restrict pscore,
  ind_t       * const restrict marker,
  val_t       * const restrict tmpmax,
  val_t       * const restrict tmpsum,
  val_t       * const restrict tmplen,
  val_t       * const restrict tmpspa,
  struct cand * const restrict tmpcnd
)
{
  ind_t ncnt = 0;

  for (ind_t j = ia[i]; j < ia[i + 1]; j++) {
    tmpmax[ja[j]] = max[j];
    tmpsum[ja[j]] = sum[j];
    tmplen[ja[j]] = len[j];
    tmpspa[ja[j]] = a[j];
  }

  for (ind_t j = 0; j < cnt; j++) {
    ind_t const k = tmpcnd[j].ind;
    val_t       s = tmpcnd[j].sim;
    ind_t const m = marker[k];
    val_t       remsim = minsim - s;

    /* reset markers to unknown value */
    marker[k] = UNKNOWN;

    /* ... */
    if (PRUNED == m) {
      /* update global counter */
      apss_nprun++;

      continue;
    }

    /* -----------------------------------------------------------------------*/
    if ( pscore[k]                           < remsim   /* Anastasiu */
      || rs1[ka[k] - 1]                      < remsim   /* Bayardo   */
      || max[ia[i + 1] - 1] * sum[ka[k] - 1] < remsim   /* Awekar    */
      || max[ka[k] - 1] * sum[ia[i + 1] - 1] < remsim ) /* ...       */
      continue;

    /* -----------------------------------------------------------------------*/
    for (ind_t jjp1 = ka[k]; jjp1 > ia[k]; jjp1--) {
      ind_t const jj = jjp1 - 1;
      if (tmpspa[ja[jj]] > 0.0) {
        remsim = minsim - s;

#if 0
        if ( tmplen[ja[jj]] * len[jj] < remsim ) /* Anastasiu */
#else
        if ( tmplen[ja[jj]] * len[jj] < remsim   /* Anastasiu */
          || rs1[jj]                  < remsim   /* Bayardo   */
          || tmpmax[ja[jj]] * sum[jj] < remsim   /* Awekar    */
          || max[jj] * tmpsum[ja[jj]] < remsim ) /* ...       */
#endif
        {
          s = 0.0;
          break;
        }

        s += tmpspa[ja[jj]] * a[jj];

        /* update global counter */
        apss_nmacs2++;
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

  for (ind_t j = ia[i]; j < ia[i + 1]; j++) {
    tmpmax[ja[j]] = 0.0;
    tmpsum[ja[j]] = 0.0;
    tmplen[ja[j]] = 0.0;
    tmpspa[ja[j]] = 0.0;
  }

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
  val_t  const * const m_max = pp->m_max;
  val_t  const * const m_sum = pp->m_sum;
  val_t  const * const m_len = pp->m_len;
  val_t  const * const m_rs1 = pp->m_rs1;
  val_t  const * const i_max = pp->i_max;
  val_t  const * const i_sum = pp->i_sum;
  val_t  const * const i_len = pp->i_len;
  val_t  const * const i_rs1 = pp->i_rs1;
  ind_t  const * const psplit = pp->psplit;
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
  val_t       * const tmpsum = GC_calloc(m_nc, sizeof(*tmpsum));
  val_t       * const tmpmax = GC_calloc(m_nc, sizeof(*tmpmax));
  val_t       * const tmplen = GC_calloc(m_nc, sizeof(*tmplen));
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
    ind_t const cnt = generate(minsim, i, m_ia, m_ja, m_ra, m_a, m_max, m_sum,
                               m_len, i_ia, i_ja, i_a, i_max, i_sum, i_len,
                               i_rs1, psplit, marker, tmpcnd);

    /* verify candidate vectors */
    ind_t const ncnt = verify(minsim, cnt, i, m_ia, m_ja, m_ka, m_a, m_max,
                              m_sum, m_len, m_rs1, pscore, marker, tmpmax,
                              tmpsum, tmplen, tmpspa, tmpcnd);

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
  GC_free(tmpsum);
  GC_free(tmpmax);
  GC_free(tmplen);
  GC_free(tmpspa);
  GC_free(tmpcnd);

  return 0;
}
