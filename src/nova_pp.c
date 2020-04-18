/* SPDX-License-Identifier: MIT */
#include <math.h>
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/nova.h"
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

    /* TODO: There is some potential here to use the Awekar bound to reduce the
     * number of items indexed --- we will need to keep track of the sum up to
     * the given point for this row and the maximum sum up to the given point
     * for all rows after this one; we will also need to keep track of the
     * maximum value up to the given point for this row and the maximum value up
     * to the given point for all rows after this one (this is different than
     * colmax).
     *
     * This will replace Bayardo bound */
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
  ind_t const nc,
  ind_t const * const ia,
  ind_t const * const ja,
  val_t const * const a,
  ind_t       * const ka,
  val_t       * const rs1,
  val_t       * const rs3,
  val_t       * const sum,
  val_t       * const max,
  val_t       * const rowmax,
  val_t       * const pfxmax,
  val_t       * const pscore
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* allocate scratch memory */
  val_t * const colmax = GC_calloc(nc, sizeof(*colmax));

  for (ind_t ip1 = nr; ip1 > 0; ip1--) {
    ind_t i = ip1 - 1;

    /* compute row maximums */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > rowmax[i])
        rowmax[i] = a[j];

    /* find prefix split for each row_i */
    ka[i] = ia[i] + filter_row(minsim, rowmax[i], ia[i + 1] - ia[i], ja + ia[i],
                               a + ia[i], colmax);

    /* XXX: (improvement) update column maximums */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > colmax[ja[j]])
        colmax[ja[j]] = a[j];

    /* compute prefix sums */
    sum[ia[i]] = a[ia[i]];
    for (ind_t j = ia[i] + 1; j < ia[i + 1]; j++)
      sum[j] = sum[j - 1] + a[j];

    /* compute prefix maximums */
    max[ia[i]] = a[ia[i]];
    for (ind_t j = ia[i] + 1; j < ia[i + 1]; j++)
      max[j] = max(max[j - 1], a[j]);

    /* XXX: Keeping column prefix max and sum across all columns could be a good
     * use for a fenwick tree */

    /* compute prefix maximums */
    for (ind_t j = ia[i]; j < ka[i]; j++)
      if (a[j] > pfxmax[i])
        pfxmax[i] = a[j];

    /* compute pscores */
    val_t b1 = 0.0, bt = 0.0, b3 = 0.0;
    for (ind_t j = ia[i]; j < ka[i]; j++) {
      rs1[j] = b1;
      rs3[j] = b3;

      b1 += a[j] * min(rowmax[i], colmax[ja[j]]);
      bt += a[j] * a[j];
      b3  = sqrtv(bt);
    }
    pscore[i] = min(b1, b3);
  }

  /* reset colmax */
  memset(colmax, 0, nc * sizeof(*colmax));

  for (ind_t i = 0; i < nr; i++) {
    val_t b1 = 0.0, bt = 1.0, b3 = 1.0;

    /* precompute b1 */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      b1 += a[j] * colmax[ja[j]];

    /* compute remscores */
    for (ind_t jp1 = ia[i + 1]; jp1 > ia[i]; jp1--) {
      ind_t const j = jp1 - 1;

      b1 -= a[j] * colmax[ja[j]];
      bt -= a[j] * a[j];
      b3  = bt < 0.0 ? 0.0 : sqrtv(bt);

      rs1[j] = b1;
      rs3[j] = b3;
    }

    /* update column maximums */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > colmax[ja[j]])
        colmax[ja[j]] = a[j];
  }

  /* free scratch memory */
  GC_free(colmax);

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
  ind_t const * const restrict ka,
  val_t const * const restrict acsr,
  val_t const * const restrict sum,
  val_t const * const restrict max,
  ind_t       * const restrict ia1,
  ind_t       * const restrict ja1,
  val_t       * const restrict acsc,
  val_t       * const restrict rs31,
  val_t       * const restrict sum1,
  val_t       * const restrict max1
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
    val_t rst = 0.0;
    for (ind_t j = ia[i]; j < ka[i]; j++)
      rst += acsr[j] * acsr[j];
    for (ind_t j = ka[i]; j < ia[i + 1]; j++) {
      ja1[ia1[ja[j]]]    = i;
      acsc[ia1[ja[j]]]   = acsr[j];
      rs31[ia1[ja[j]]]   = sqrtv(rst);
      sum1[ia1[ja[j]]]   = sum[j];
      max1[ia1[ja[j]]++] = max[j];
      rst += acsr[j] * acsr[j];
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
  val_t  const * const m_sum,
  val_t  const * const m_max,
  Matrix       * const I,
  val_t ** const i_rs3,
  val_t ** const i_sum,
  val_t ** const i_max
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* unpack /M/ */
  ind_t const m_nr  = M->nr;
  ind_t const m_nc  = M->nc;
  ind_t const * const m_ia  = M->ia;
  ind_t const * const m_ja  = M->ja;
  val_t const * const m_a   = M->a;

  /* compute the size of the inverted index */
  ind_t i_nnz = 0;
  for (ind_t i = 0; i < m_nr; i++)
    i_nnz += m_ia[i + 1] - m_ka[i];

  /* allocate memory for inverted index */
  ind_t * const i_ia  = GC_malloc((m_nc + 1) * sizeof(*i_ia));
  ind_t * const i_ja  = GC_malloc(i_nnz * sizeof(*i_ja));
  val_t * const i_a   = GC_malloc(i_nnz * sizeof(*i_a));
               *i_rs3 = GC_malloc(i_nnz * sizeof(*i_rs3));
               *i_sum = GC_malloc(i_nnz * sizeof(*i_sum));
               *i_max = GC_malloc(i_nnz * sizeof(*i_max));

  csrcsc(m_nr, m_nc, m_ia, m_ja, m_ka, m_a, m_sum, m_max, i_ia, i_ja, i_a,
         *i_rs3, *i_sum, *i_max);

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
/*! Precompute row starts for inverted index. */
/*----------------------------------------------------------------------------*/
static int
diidx(
  val_t  const         minsim,
  Matrix const * const M,
  ind_t        * const m_ra,
  Matrix const * const I,
  val_t  const * const rowmax
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* unpack /M/ */
  ind_t const m_nr  = M->nr;
  ind_t const m_nc  = M->nc;
  ind_t const * const m_ia = M->ia;
  ind_t const * const m_ja = M->ja;

  /* unpack /I/ */
  ind_t const * const i_ia = I->ia;
  ind_t const * const i_ja = I->ja;

  /* allocate scratch memory */
  ind_t * const i_ra = GC_malloc(m_nc * sizeof(*i_ra));

  /* initialize /i_ra/ */
  memcpy(i_ra, i_ia, m_nc * sizeof(*i_ra));

  /* remove rows from index that are too short */
  for (ind_t i = 0; i < m_nr; i++) {
    val_t const sz1 = minsim / rowmax[i];

    for (ind_t j = m_ia[i]; j < m_ia[i + 1]; j++) {
      ind_t const jj = m_ja[j];

      for (; i_ra[jj] < i_ia[jj + 1]; i_ra[jj]++) {
        ind_t const k = i_ja[i_ra[jj]];
        if ((m_ia[k + 1] - m_ia[k]) * rowmax[k] >= sz1)
          break;
      }

      m_ra[j] = i_ra[jj];
    }
  }

  /* free scratch memory */
  GC_free(i_ra);

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
  free(pp->ra);
  free(pp->m_rs1);
  free(pp->m_rs3);
  free(pp->m_sum);
  free(pp->m_max);
  free(pp->i_rs3);
  free(pp->i_sum);
  free(pp->i_max);
  free(pp->rowmax);
  free(pp->pfxmax);
  free(pp->pscore);
}

/*----------------------------------------------------------------------------*/
/*! Function to preprocess matrix for APSS. */
/*----------------------------------------------------------------------------*/
EFIKA_APSS_EXPORT int
apss_nova_pp(
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
  ind_t const m_nr  = M->nr;
  ind_t const m_nc  = M->nc;
  ind_t const m_nnz = M->nnz;
  ind_t const * const m_ia = M->ia;
  ind_t const * const m_ja = M->ja;
  val_t const * const m_a  = M->a;

  /* unpack /pp/ */
  Matrix * const I = &(pp->I);

  /* allocate memory for other preprocessed data */
  ind_t * const m_ka  = GC_malloc(m_nr * sizeof(*m_ka));
  ind_t * const m_ra  = GC_malloc(m_nnz * sizeof(*m_ra));
  val_t * const m_rs1 = GC_malloc(m_nnz * sizeof(*m_rs1));
  val_t * const m_rs3 = GC_malloc(m_nnz * sizeof(*m_rs3));
  val_t * const m_sum = GC_calloc(m_nnz, sizeof(*m_sum));
  val_t * const m_max = GC_calloc(m_nnz, sizeof(*m_max));
  val_t * const rowmax = GC_calloc(m_nr, sizeof(*rowmax));
  val_t * const pfxmax = GC_calloc(m_nr, sizeof(*pfxmax));
  val_t * const pscore = GC_calloc(m_nr, sizeof(*pscore));
  val_t * i_rs3;
  val_t * i_sum;
  val_t * i_max;

  /* init /I/ */
  err = Matrix_init(I);
  GC_assert(!err);

  /* precompute row prefixes and their statistics */
  err = rfilt(minsim, m_nr, m_nc, m_ia, m_ja, m_a, m_ka, m_rs1, m_rs3, m_sum,
              m_max, rowmax, pfxmax, pscore);
  GC_assert(!err);

  /* precompute dynamic inverted index */
  err = fiidx(M, m_ka, m_sum, m_max, I, &i_rs3, &i_sum, &i_max);
  GC_assert(!err);

  /* precompute dynamic inverted index */
  err = diidx(minsim, M, m_ra, I, rowmax);
  GC_assert(!err);

  /* record info in /pp/ */
  pp->ka    = m_ka;
  pp->ra    = m_ra;
  pp->m_rs1 = m_rs1;
  pp->m_rs3 = m_rs3;
  pp->m_sum = m_sum;
  pp->m_max = m_max;
  pp->i_rs3 = i_rs3;
  pp->i_sum = i_sum;
  pp->i_max = i_max;
  pp->rowmax = rowmax;
  pp->pfxmax = pfxmax;
  pp->pscore = pscore;

  /* record payload in /M/ */
  M->pp = pp;
  M->pp_free = &pp_free;

  return 0;
}
