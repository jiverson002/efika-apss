/* SPDX-License-Identifier: MIT */
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
/* Precompute prefix for a single row. */
/*----------------------------------------------------------------------------*/
static inline ind_t
filter_row(
  val_t const minsim,
  ind_t const n,
  ind_t const * const ja,
  val_t const * const a,
  val_t const * const tmpmax
)
{
  val_t b1 = 0.0, mx = 0.0;

  /* compute row maximum */
  for (ind_t j = 0; j < n; j++)
    if (a[j] > mx)
      mx = a[j];

  for (ind_t j = 0; j < n; j++) {
    /* Bayardo bound */
    b1 += a[j] * min(mx, tmpmax[ja[j]]);

    if (b1 >= minsim) {
      /* return prefix split */
      return j;
    }
  }

  /* return end-of-row */
  return n;
}

/*----------------------------------------------------------------------------*/
/*! Precompute row prefixes. */
/*----------------------------------------------------------------------------*/
static int
filter(
  val_t const minsim,
  ind_t const nr,
  ind_t const nc,
  ind_t const * const ia,
  ind_t const * const ja,
  val_t const * const a,
  ind_t       * const ka
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* allocate scratch memory */
  val_t * const tmpmax = GC_calloc(nc, sizeof(*tmpmax));

  /* compute column maximums */
  for (ind_t i = 0; i < nr; i++)
    for(ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > tmpmax[ja[j]])
        tmpmax[ja[j]] = a[j];

  /* find the prefixes in the matrix */
  for (ind_t i = 0; i < nr; i++)
    /* find prefix split for each row_i */
    ka[i] = ia[i] + filter_row(minsim, ia[i + 1] - ia[i], ja + ia[i], a + ia[i],
                               tmpmax);

  /* free scratch memory */
  GC_free(tmpmax);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*! Converts a sparse matrix in the CSR format to the CSC format. */
/*----------------------------------------------------------------------------*/
static void
csrcsc(
  ind_t const nr,
  ind_t const nc,
  ind_t const * const restrict ia,
  ind_t const * const restrict ja,
  val_t const * const restrict acsr,
  ind_t const * const restrict ka,
  ind_t       * const restrict ia1,
  ind_t       * const restrict ja1,
  val_t       * const restrict acsc
)
{
  memset(ia1, 0, (nc + 1) * sizeof(*ia1));

  for (ind_t i = 0; i < nr; i++)
    for (ind_t j = ka[i]; j < ia[i + 1]; j++)
      ia1[ja[j]]++;

  for (ind_t i = 0, p = 0; i <= nc; i++) {
    ind_t const t = ia1[i];
    ia1[i] = p;
    p += t;
  }

  for (ind_t i = 0; i < nr; i++) {
    for (ind_t j = ka[i]; j < ia[i + 1]; j++) {
      ja1[ia1[ja[j]]]    = i;
      acsc[ia1[ja[j]]++] = acsr[j];
    }
  }

  for (ind_t i = nc; i > 0; i--)
    ia1[i] = ia1[i - 1];
  ia1[0] = 0;
}

/*----------------------------------------------------------------------------*/
/*! Create inverted index based on precomputed row prefixes. */
/*----------------------------------------------------------------------------*/
static int
fiidx(
  val_t  const         minsim,
  Matrix const * const M,
  ind_t        * const m_ka,
  Matrix       * const I
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* unpack /M/ */
  ind_t const         m_nr  = M->nr;
  ind_t const         m_nc  = M->nc;
  ind_t const * const m_ia  = M->ia;
  ind_t const * const m_ja  = M->ja;
  val_t const * const m_a   = M->a;

  /* precompute row prefixes */
  int err = filter(minsim, m_nr, m_nc, m_ia, m_ja, m_a, m_ka);
  GC_assert(!err);

  /* compute the size of the inverted index */
  ind_t i_nnz = 0;
  for (ind_t i = 0; i < m_nr; i++)
    i_nnz += m_ia[i + 1] - m_ka[i];

  /* allocate memory for inverted index */
  ind_t * const i_ia = GC_malloc((m_nc + 1) * sizeof(*i_ia));
  ind_t * const i_ja = GC_malloc(i_nnz * sizeof(*i_ja));
  val_t * const i_a  = GC_malloc(i_nnz * sizeof(*i_a));

  csrcsc(m_nr, m_nc, m_ia, m_ja, m_a, m_ka, i_ia, i_ja, i_a);

  /* record relevant info in /I/ */
  I->nr  = m_nc;
  I->nc  = m_nr;
  I->nnz = i_nnz;
  I->ia  = i_ia;
  I->ja  = i_ja;
  I->a   = i_a;

  return 0;
}

/*----------------------------------------------------------------------------*/
/*! Function to free any memory allocated during preprocessing. */
/*----------------------------------------------------------------------------*/
static void
pp_free(
  void * const ptr
)
{
  struct pp_payld * pp = (struct pp_payld *)ptr;
  Matrix_free(&(pp->I));
  free(pp->ka);
}

/*----------------------------------------------------------------------------*/
/*! Function to preprocess matrix for APSS. */
/*----------------------------------------------------------------------------*/
EFIKA_APSS_EXPORT int
apss_allpairs_pp(
  val_t const minsim,
  Matrix * const M
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* reorder columns in decreasing order of degree */
  int err = Matrix_cord(M, DEG | DSC);
  GC_assert(!err);

  /* reorder rows in decreasing order of row maximum */
  err = Matrix_rord(M, VAL | DSC);
  GC_assert(!err);

  /* reorder each row in increasing order of column id */
  err = Matrix_sort(M, COL | ASC);
  GC_assert(!err);

  /* ... */
  struct pp_payld * pp = GC_malloc(sizeof(*pp));

  /* unpack /M/ */
  ind_t const nr = M->nr;

  /* unpack /pp/ */
  Matrix * const I = &(pp->I);

  /* init /ka/ */
  ind_t * const ka = GC_malloc(nr * sizeof(*ka));

  /* init /I/ */
  err = Matrix_init(I);
  GC_assert(!err);

  /* ... */
  err = fiidx(minsim, M, ka, I);
  GC_assert(!err);

  /* record info in /pp/ */
  pp->ka = ka;

  /* record payload in /M/ */
  M->pp = pp;
  M->pp_free = &pp_free;

  return 0;
}
