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

#include <stdio.h>

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
  val_t const * const restrict m_len,
  ind_t const * const restrict i_ia,
  ind_t const * const restrict i_ja,
  val_t const * const restrict i_a,
  val_t const * const restrict i_len,
  ind_t       * const restrict marker,
  struct cand * const restrict tmpcnd
)
{
  ind_t cnt = 0;
  val_t rs1 = 1.0;

  /* iterate through each non-zero column, j, for this row */
  for (ind_t jj = m_ia[i + 1]; jj > m_ia[i]; jj--) {
    ind_t const j    = m_ja[jj - 1];
    val_t const v    = m_a[jj - 1];
    ind_t const i_ra = m_ra[jj - 1];

    /* retrieve precomputed values */
    val_t const len = m_len[jj - 1];

    /* determine if any new candidates should be allowed */
    bool const allow_unknown  = /* jj > m_ka[i] && */ min(rs1, len) >= minsim;

    /* retrieve precomputed values */
    rs1 = m_rs1[jj - 1];

    /* iterate through the rows, k,  that have indexed column j */
    for (ind_t kk = i_ra; kk < i_ia[j + 1]; kk++) {
      ind_t const k = i_ja[kk];
      val_t const w = i_a[kk];
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

          /* update global counter */
          apss_nmacs1++;
        }
        break;

        case PRUNED:
        break;

        default:
        /* Anastasiu pruning */
        if (tmpcnd[m].sim + len * i_len[kk] < minsim) {
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

  (void)m_ka;
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
  val_t const * const restrict sum,
  val_t const * const restrict max,
  val_t const * const restrict len,
  val_t const * const restrict pscore,
  ind_t       * const restrict marker,
  val_t       * const restrict tmpsum,
  val_t       * const restrict tmpmax,
  val_t       * const restrict tmplen,
  val_t       * const restrict tmpspa,
  struct cand * const restrict tmpcnd
)
{
  ind_t ncnt = 0;
  ind_t const ln = ia[i + 1] - ia[i];
  val_t const mx = max[ia[i + 1] - 1];

  BLAS_vsctr(ia[i + 1] - ia[i],   a + ia[i], ja + ia[i], tmpspa);
  BLAS_vsctr(ia[i + 1] - ia[i], sum + ia[i], ja + ia[i], tmpsum);
  BLAS_vsctr(ia[i + 1] - ia[i], max + ia[i], ja + ia[i], tmpmax);
  BLAS_vsctr(ia[i + 1] - ia[i], len + ia[i], ja + ia[i], tmplen);

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

    /* -----------------------------------------------------------------------*/
    val_t const ub1 = min(
      min(ln, ka[k] - ia[k]) * mx * max[ka[k] - 1], /* Bayardo */
      pscore[k]                                     /* pscore  */
    );
    if (s + ub1 < minsim)
      continue;

#if 0
    /* -----------------------------------------------------------------------*/
    /* fast-forward to next matching column */
    ind_t jjp1;
    for (jjp1 = ka[k]; jjp1 > ia[k] && tmpspa[ja[jjp1 - 1]] == 0.0; jjp1--);
    ind_t const jj = jjp1 - 1;

    /* no more matching columns */
    if (jjp1 == ia[k] && s < minsim)
      continue;

    val_t const ub2 = min(
      tmpmax[ja[jj]] * sum[jj], /* Awekar    */
      max[jj] * tmpsum[ja[jj]]  /* ...       */
      tmplen[ja[jj]] + len[jj]  /* Anastasiu */
    );
    if (s + ub2 < minsim)
      continue;
#endif

    /* -----------------------------------------------------------------------*/
    /* Anastasiu dot product */
    for (ind_t jjp1 = ka[k]; jjp1 > ia[k]; jjp1--) {
      ind_t const jj = jjp1 - 1;
      if (tmpspa[ja[jj]] > 0.0) {
#if 1
        val_t const ub3 = min3(
          tmpmax[ja[jj]] * sum[jj], /* Awekar    */
          max[jj] * tmpsum[ja[jj]], /* ...       */
          tmplen[ja[jj]] + len[jj]  /* Anastasiu */
        );
        if (s + ub3 < minsim) {
          s = 0.0;
          break;
        }
#endif

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

  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmpsum);
  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmpmax);
  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmplen);
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
  val_t  const * const m_sum = pp->m_sum;
  val_t  const * const m_max = pp->m_max;
  val_t  const * const m_len = pp->m_len;
  val_t  const * const i_len = pp->i_len;
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
    ind_t const cnt = generate(minsim, i, m_ia, m_ja, m_ka, m_ra, m_a, m_rs1,
                               m_len, i_ia, i_ja, i_a, i_len, marker, tmpcnd);

    /* verify candidate vectors */
    ind_t const ncnt = verify(minsim, cnt, i, m_ia, m_ja, m_ka, m_a, m_sum,
                              m_max, m_len, pscore, marker, tmpsum, tmpmax,
                              tmplen, tmpspa, tmpcnd);

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
