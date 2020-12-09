/* SPDX-License-Identifier: MIT */
#include <stdbool.h>
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/mmjoin.h"
#include "efika/apss/rename.h"
#include "efika/core/blas.h"
#include "efika/core/export.h"
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
  ind_t const * const restrict ia,
  ind_t const * const restrict ja,
  val_t const * const restrict a,
  val_t       * const restrict l,
  ind_t const * const restrict ia1,
  ind_t const * const restrict ja1,
  ind_t       * const restrict ra1,
  val_t const * const restrict a1,
  val_t const * const restrict l1,
  val_t const * const restrict rowmax,
  val_t const * const restrict colmax,
  ind_t       * const restrict marker,
  val_t       * const restrict tmpl,
  struct cand * const restrict tmpcnd
)
{
  ind_t cnt = 0;
  val_t rs1 = 0.0, rs2 = 1.0, lx = 0.5;

  /* precompute the minimum size */
  val_t const sz1 = minsim / rowmax[i];

  /* precompute maximum dot product */
  for (ind_t jj = ia[i]; jj < ia[i + 1]; jj++)
    rs1 += a[jj] * colmax[ja[jj]];

  /* iterate through each non-zero column, j, for this row */
  for (ind_t jj = ia[i + 1]; jj > ia[i]; jj--) {
    ind_t const j = ja[jj - 1];
    val_t const v = a[jj - 1];

    /* ... */
    bool const allow_unknown = min(rs1, rs2) >= minsim;

    /* remove rows from index that are too short */
    for (; ra1[j] < ia1[j + 1]; ra1[j]++)
      if ((ia[ja1[ra1[j]] + 1] - ia[ja1[ra1[j]]]) * rowmax[ja1[ra1[j]]] >= sz1)
        break;

    /* ... */
    lx -= 0.5 * v * v;

    /* iterate through the rows, k,  that have indexed column j */
    for (ind_t kk = ra1[j]; kk < ia1[j + 1]; kk++) {
      ind_t const k  = ja1[kk];
      val_t const w  = a1[kk];
      val_t const ly = l1[kk];
      ind_t       m  = marker[k];

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

        /* populate marker */
        m = marker[k] = cnt++;

        /* fall-through */

        default:
        /* update partial dot product for candidate row */
        tmpcnd[m].sim += v * w;

        /* update global counter */
        apss_nmacs1++;

        /* Lee pruning */
        if (tmpcnd[m].sim + lx + ly < minsim)
          marker[k] = PRUNED;
      }
    }

    /* ... */
    rs1 -= v * colmax[j];
    rs2 -= 0.5 * v * v;
    l[jj - 1] = lx;
    tmpl[j] = lx;
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
  val_t const * const restrict l,
  val_t const * const restrict rowmax,
  val_t const * const restrict pfxmax,
  ind_t       * const restrict marker,
  val_t       * const restrict tmpl,
  val_t       * const restrict tmpspa,
  struct cand * const restrict tmpcnd
)
{
  ind_t ncnt = 0, len = ia[i + 1] - ia[i];
  val_t const mx = rowmax[i];

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
    if (s + min(ka[k] - ia[k], len) * mx * pfxmax[k] < minsim)
      continue;

    /* Lee dot product */
    ind_t unmatched = 0;
    for(ind_t jj = ka[k]; jj > ia[k]; jj--){
      ind_t const jjj = jj - 1;
      if (tmpspa[ja[jjj]] > 0.0) {
        s += tmpspa[ja[jjj]] * a[jjj];

        /* update global counter */
        apss_nmacs2++;

        if (unmatched > 1) {
          if (s + tmpl[ja[jjj]] + l[jjj] < minsim) {
            s = 0.0;
            break;
          }
          unmatched = 0;
        }
      } else {
        unmatched++;
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

  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmpl);
  BLAS_vsctrz(ia[i + 1] - ia[i], ja + ia[i], tmpspa);

  return ncnt;
}

/*----------------------------------------------------------------------------*/
/*! Function to compute APSS of a matrix. */
/*----------------------------------------------------------------------------*/
EFIKA_EXPORT int
apss_mmjoin(
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
    int err = apss_mmjoin_pp(minsim, M);
    GC_assert(!err);
  }

  /* unpack /M/ */
  ind_t const m_nr  = M->nr;
  ind_t const m_nc  = M->nc;
  ind_t const m_nnz = M->nnz;
  ind_t const * const m_ia = M->ia;
  ind_t const * const m_ja = M->ja;
  val_t const * const m_a  = M->a;
  struct pp_payld const * const pp = M->pp;

  /* unpack /pp/ */
  Matrix const * const I    = &(pp->I);
  ind_t  const * const m_ka = pp->ka;
  val_t  const * const i_l  = pp->l;
  val_t  const * const rowmax = pp->rowmax;
  val_t  const * const colmax = pp->colmax;
  val_t  const * const pfxmax = pp->pfxmax;

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
  val_t       * const m_l    = GC_malloc(m_nnz * sizeof(*m_l));
  ind_t       * const i_ra   = GC_malloc(m_nc * sizeof(*i_ra));
  ind_t       * const marker = GC_malloc(m_nr * sizeof(*marker));
  val_t       * const tmpl   = GC_calloc(m_nc, sizeof(*tmpl));
  val_t       * const tmpspa = GC_calloc(m_nc, sizeof(*tmpspa));
  struct cand * const tmpcnd = GC_malloc(m_nr * sizeof(*tmpcnd));

  /* make a copy of preprocessed data */
  memcpy(i_ra, i_ia, m_nc * sizeof(*i_ra));

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
    ind_t const cnt = generate(minsim, i, m_ia, m_ja, m_a, m_l, i_ia, i_ja,
                               i_ra, i_a, i_l, rowmax, colmax, marker, tmpl,
                               tmpcnd);

    /* verify candidate vectors */
    ind_t const ncnt = verify(minsim, cnt, i, m_ia, m_ja, m_ka, m_a, m_l,
                              rowmax, pfxmax, marker, tmpl, tmpspa, tmpcnd);

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
    free(M->pp);
    M->pp = NULL;
  }

  /* free scratch memory */
  GC_free(m_l);
  GC_free(i_ra);
  GC_free(marker);
  GC_free(tmpl);
  GC_free(tmpspa);
  GC_free(tmpcnd);

  return 0;
}
