/* SPDX-License-Identifier: MIT */
#include <math.h>
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/l2ap.h"
#include "efika/apss/export.h"
#include "efika/apss/rename.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"

/*----------------------------------------------------------------------------*/
/* Precompute prefix for a single row. */
/*----------------------------------------------------------------------------*/
static inline ind_t
filter_row(
  val_t const minsim,
  val_t const mx,
  ind_t const n,
  ind_t const * const ja,
  val_t const * const a,
  val_t const * const colmax
)
{
  ind_t j;
  val_t b1 = 0.0, bt = 0.0, b3 = 0.0;

  for (j = 0; j < n && min(b1, b3) < minsim; j++) {
    /* Bayardo bound */
    b1 += a[j] * min(mx, colmax[ja[j]]);

   /* Anastasiu bound */
    bt += a[j] * a[j];
    b3  = sqrtv(bt);

    //pscore[ii] = min(b1, b3);
  }

  /* return prefix split */
  return j - 1;
}

/*----------------------------------------------------------------------------*/
/*! Precompute row prefixes. */
/*----------------------------------------------------------------------------*/
static int
rfilt(
  val_t const minsim,
  ind_t const nr,
  ind_t const * const ia,
  ind_t const * const ja,
  val_t const * const a,
  ind_t       * const ka,
  val_t       * const rowmax,
  val_t       * const colmax,
  val_t       * const pfxmax,
  val_t       * const pscore
)
{
  /* compute column maximums */
  for (ind_t i = 0; i < nr; i++)
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > colmax[ja[j]])
        colmax[ja[j]] = a[j];

  /* find the prefixes in the matrix */
  for (ind_t i = 0; i < nr; i++) {
    /* compute row maximums */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > rowmax[i])
        rowmax[i] = a[j];

    /* find prefix split for each row_i */
    ka[i] = ia[i] + filter_row(minsim, rowmax[i], ia[i + 1] - ia[i], ja + ia[i],
                               a + ia[i], colmax);

    /* compute prefix maximums */
    for (ind_t j = ia[i]; j < ka[i]; j++)
      if (a[j] > pfxmax[i])
        pfxmax[i] = a[j];

    /* compute pscores */
    val_t b1 = 0.0, bt = 0.0;
    for (ind_t j = ia[i]; j < ka[i]; j++) {
      b1 += a[j] * min(rowmax[i], colmax[ja[j]]);
      bt += a[j] * a[j];
    }
    pscore[i] = min(b1, sqrtv(bt));
  }

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
  val_t       * const restrict acsc,
  val_t       * const restrict l1
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
    val_t l = 0.0;
    for (ind_t j = ia[i]; j < ka[i]; j++)
      l += acsr[j] * acsr[j];
    for (ind_t j = ka[i]; j < ia[i + 1]; j++) {
      ja1[ia1[ja[j]]]  = i;
      acsc[ia1[ja[j]]] = acsr[j];
      l1[ia1[ja[j]]++] = sqrtv(l);
      l += acsr[j] * acsr[j];
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
  Matrix const * const M,
  ind_t  const * const m_ka,
  Matrix       * const I,
  val_t ** const i_l
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

  /* compute the size of the inverted index */
  ind_t i_nnz = 0;
  for (ind_t i = 0; i < m_nr; i++)
    i_nnz += m_ia[i + 1] - m_ka[i];

  /* allocate memory for inverted index */
  ind_t * const i_ia = GC_malloc((m_nc + 1) * sizeof(*i_ia));
  ind_t * const i_ja = GC_malloc(i_nnz * sizeof(*i_ja));
  val_t * const i_a  = GC_malloc(i_nnz * sizeof(*i_a));
               *i_l  = GC_malloc(i_nnz * sizeof(*i_l));

  csrcsc(m_nr, m_nc, m_ia, m_ja, m_a, m_ka, i_ia, i_ja, i_a, *i_l);

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
  free(pp->l);
  free(pp->rowmax);
  free(pp->pfxmax);
  free(pp->pscore);
}

/*----------------------------------------------------------------------------*/
/*! Function to preprocess matrix for APSS. */
/*----------------------------------------------------------------------------*/
EFIKA_APSS_EXPORT int
apss_l2ap_pp(
  val_t const minsim,
  Matrix * const M
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* ... */
  if (!M)
    return -1;

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
  ind_t const m_nr = M->nr;
  ind_t const m_nc = M->nc;
  ind_t const * const m_ia  = M->ia;
  ind_t const * const m_ja  = M->ja;
  val_t const * const m_a   = M->a;

  /* unpack /pp/ */
  Matrix * const I = &(pp->I);

  /* allocate memory for other preprocessed data */
  ind_t * const m_ka = GC_malloc(m_nr * sizeof(*m_ka));
  val_t * const rowmax = GC_calloc(m_nr, sizeof(*rowmax));
  val_t * const pfxmax = GC_calloc(m_nr, sizeof(*pfxmax));
  val_t * const pscore = GC_calloc(m_nr, sizeof(*pscore));
  val_t * i_l;

  /* allocate scratch memory */
  val_t * const colmax = GC_calloc(m_nc, sizeof(*colmax));

  /* init /I/ */
  err = Matrix_init(I);
  GC_assert(!err);

  /* precompute row prefixes and their statistics */
  err = rfilt(minsim, m_nr, m_ia, m_ja, m_a, m_ka, rowmax, colmax, pfxmax,
              pscore);
  GC_assert(!err);

  /* precompute dynamic inverted index */
  err = fiidx(M, m_ka, I, &i_l);
  GC_assert(!err);

  /* record info in /pp/ */
  pp->ka = m_ka;
  pp->l  = i_l;
  pp->rowmax = rowmax;
  pp->pfxmax = pfxmax;
  pp->pscore = pscore;

  /* record payload in /M/ */
  M->pp = pp;
  M->pp_free = &pp_free;

  /* free scratch memory */
  GC_free(colmax);

  return 0;
}
