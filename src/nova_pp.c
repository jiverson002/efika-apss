/* SPDX-License-Identifier: MIT */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/bit.h"
#include "efika/apss/nova.h"
#include "efika/apss/export.h"
#include "efika/apss/rename.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"

/*----------------------------------------------------------------------------*/
/* Row prefix statistics. */
/*----------------------------------------------------------------------------*/
typedef struct {
  ind_t n;
  val_t pscore;
} Stat;

/*----------------------------------------------------------------------------*/
/* Precompute prefix for a single row. */
/*----------------------------------------------------------------------------*/
static inline Stat
filt_row(
  val_t const minsim,
  ind_t const n,
  ind_t const * const ja,
  val_t const * const a,
  val_t       * const rowmax,
  val_t       * const rowsum,
  val_t       * const rowsqr,
  val_t       * const rowlen,
  val_t const * const colmax,
  val_t const * const colmax_,
  val_t const * const colsum_,
  val_t const * const colsqr_,
  val_t const * const collen_
)
{
  rowmax[0] = a[0];        /* Awekar    */
  rowsum[0] = a[0];        /* ...       */
  rowsqr[0] = a[0] * a[0]; /* Anastasiu */
  rowlen[0] = a[0];        /* ...       */

  /* XXX: Due to rounding error, this could actually be slightly above 1.0.
   *      This could easily be corrected with:
   *
   *      clen = clen <= 1.0 ? clen : 1.0;
   */
  val_t clen = BIT_sum(ja[0], collen_, BIT_op_max);

  /* Bayardo bound */
  val_t b1 = a[0] * colmax[ja[0]];

  /* dot-product upper-bound */
  val_t ub = min(
    b1,                      /* Bayardo  */
    rowlen[0] * clen         /* Anastasiu */
  );

  ind_t j = 1;
  val_t pscore = 0.0;
  for (; j < n && ub < minsim; j++) {
    rowmax[j] = BIT_op_max(rowmax[j - 1], a[j]);
    rowsum[j] = BIT_op_sum(rowsum[j - 1], a[j]);
    rowsqr[j] = BIT_op_sum(rowsqr[j - 1], a[j] * a[j]);
    rowlen[j] = sqrtv(rowsqr[j]);

    /* record pscore before updating ub */
    pscore = ub;

    clen = BIT_sum(ja[j], collen_, BIT_op_max);

    b1 += a[j] * colmax[ja[j]];

    ub = min(
      b1,                      /* Bayardo   */
      rowlen[j] * clen         /* Anastasiu */
    );
  }

  for (ind_t jj = j; jj < n; jj++) {
    rowmax[jj] = BIT_op_max(rowmax[jj - 1], a[jj]);
    rowsum[jj] = BIT_op_sum(rowsum[jj - 1], a[jj]);
    rowsqr[jj] = BIT_op_sum(rowsqr[jj - 1], a[jj] * a[jj]);
    rowlen[jj] = sqrtv(rowsqr[jj]);
  }

  /* return prefix split */
  return (Stat){ j - 1, pscore };

  (void)colmax_;
  (void)colsum_;
  (void)colsqr_;
}

/*----------------------------------------------------------------------------*/
/* Precompute prefix statistics for a single row. */
/*----------------------------------------------------------------------------*/
static inline void
stat_row(
  ind_t const nc,
  ind_t const n,
  ind_t const * const ja,
  val_t const * const rowmax,
  val_t const * const rowsum,
  val_t const * const rowsqr,
  val_t const * const rowlen,
  val_t       * const colmax,
  val_t       * const colsum,
  val_t       * const colsqr,
  val_t       * const collen
)
{
  for (ind_t j = 0; j < n; j++) {
    BIT_add(rowmax[j], nc, ja[j], colmax, BIT_op_max);
    BIT_add(rowsum[j], nc, ja[j], colsum, BIT_op_max);
    BIT_add(rowsqr[j], nc, ja[j], colsqr, BIT_op_max);
    BIT_add(rowlen[j], nc, ja[j], collen, BIT_op_max);
  }
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
  val_t       * const max,
  val_t       * const sum,
  val_t       * const sqr,
  val_t       * const len,
  val_t       * const pscore
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* allocate scratch memory */
  val_t * const colmax  = GC_calloc(nc, sizeof(*colmax));
  val_t * const colmax_ = GC_calloc((nc + 1), sizeof(*colmax_));
  val_t * const colsum_ = GC_calloc((nc + 1), sizeof(*colsum_));
  val_t * const colsqr_ = GC_calloc((nc + 1), sizeof(*colsqr_));
  val_t * const collen_ = GC_calloc((nc + 1), sizeof(*collen_));

  for (ind_t ip1 = nr; ip1 > 0; ip1--) {
    ind_t i = ip1 - 1;

    /* find prefix split and pscore for each row_i */
    Stat const pfx = filt_row(minsim, ia[i + 1] - ia[i], ja + ia[i], a + ia[i],
                              max + ia[i], sum + ia[i], sqr + ia[i], len +
                              ia[i], colmax, colmax_, colsum_, colsqr_,
                              collen_);
    ka[i] = ia[i] + pfx.n;
    pscore[i] = pfx.pscore;

    /* update global stats for each row_i */
    stat_row(nc, ia[i + 1] - ia[i], ja + ia[i], max + ia[i], sum + ia[i],
             sqr + ia[i], len + ia[i], colmax_, colsum_, colsqr_, collen_);

    /* update column maximums */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > colmax[ja[j]])
        colmax[ja[j]] = a[j];
  }

  /* reset colmax */
  memset(colmax, 0, nc * sizeof(*colmax));

  for (ind_t i = 0; i < nr; i++) {
    /* compute final remscore */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      rs1[ia[i + 1] - 1] += a[j] * colmax[ja[j]];

    /* compute remscores */
    for (ind_t jp1 = ia[i + 1] - 1; jp1 > ia[i]; jp1--)
      rs1[jp1 - 1] = rs1[jp1] - a[jp1] * colmax[ja[jp1]];

    /* update column maximums */
    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      if (a[j] > colmax[ja[j]])
        colmax[ja[j]] = a[j];
  }

  /* free scratch memory */
  GC_free(colmax);
  GC_free(colmax_);
  GC_free(colsum_);
  GC_free(colsqr_);
  GC_free(collen_);

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
  val_t const * const restrict max,
  val_t const * const restrict sum,
  val_t const * const restrict sqr,
  val_t const * const restrict len,
  ind_t       * const restrict ia1,
  ind_t       * const restrict ja1,
  val_t       * const restrict acsc,
  val_t       * const restrict max1,
  val_t       * const restrict sum1,
  val_t       * const restrict sqr1,
  val_t       * const restrict len1
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
      acsc[ia1[ja[j]]]   = acsr[j];
      max1[ia1[ja[j]]]   = max[j];
      sum1[ia1[ja[j]]]   = sum[j];
      sqr1[ia1[ja[j]]]   = sqr[j];
      len1[ia1[ja[j]]++] = len[j];
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
  val_t  const * const m_max,
  val_t  const * const m_sum,
  val_t  const * const m_sqr,
  val_t  const * const m_len,
  Matrix       * const I,
  val_t ** const i_max,
  val_t ** const i_sum,
  val_t ** const i_sqr,
  val_t ** const i_len
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
               *i_max = GC_malloc(i_nnz * sizeof(*i_max));
               *i_sum = GC_malloc(i_nnz * sizeof(*i_sum));
               *i_sqr = GC_malloc(i_nnz * sizeof(*i_sqr));
               *i_len = GC_malloc(i_nnz * sizeof(*i_len));

  csrcsc(m_nr, m_nc, m_ia, m_ja, m_ka, m_a, m_max, m_sum, m_sqr, m_len, i_ia,
         i_ja, i_a, *i_max, *i_sum, *i_sqr, *i_len);

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
  val_t  const * const m_max,
  Matrix const * const I
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
    val_t const sz1 = minsim / m_max[m_ia[i + 1] - 1];

    for (ind_t j = m_ia[i]; j < m_ia[i + 1]; j++) {
      ind_t const jj = m_ja[j];

      for (; i_ra[jj] < i_ia[jj + 1]; i_ra[jj]++) {
        ind_t const k = i_ja[i_ra[jj]];
        if ((m_ia[k + 1] - m_ia[k]) * m_max[m_ia[k + 1] - 1] >= sz1)
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

  free(pp->m_max);
  free(pp->m_sum);
  free(pp->m_sqr);
  free(pp->m_len);

  free(pp->i_max);
  free(pp->i_sum);
  free(pp->i_sqr);
  free(pp->i_len);

  free(pp->m_rs1);
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
  val_t * const m_rs1 = GC_calloc(m_nnz, sizeof(*m_rs1));
  val_t * const m_max = GC_calloc(m_nnz, sizeof(*m_max));
  val_t * const m_sum = GC_calloc(m_nnz, sizeof(*m_sum));
  val_t * const m_sqr = GC_calloc(m_nnz, sizeof(*m_sqr));
  val_t * const m_len = GC_calloc(m_nnz, sizeof(*m_len));
  val_t * const pscore = GC_calloc(m_nr, sizeof(*pscore));
  val_t * i_max;
  val_t * i_sum;
  val_t * i_sqr;
  val_t * i_len;

  /* init /I/ */
  err = Matrix_init(I);
  GC_assert(!err);

  /* precompute row prefixes and their statistics */
  err = rfilt(minsim, m_nr, m_nc, m_ia, m_ja, m_a, m_ka, m_rs1, m_max, m_sum,
              m_sqr, m_len, pscore);
  GC_assert(!err);

  /* precompute dynamic inverted index */
  err = fiidx(M, m_ka, m_max, m_sum, m_sqr, m_len, I, &i_max, &i_sum, &i_sqr,
              &i_len);
  GC_assert(!err);

  /* precompute dynamic inverted index */
  err = diidx(minsim, M, m_ra, m_max, I);
  GC_assert(!err);

  /* record info in /pp/ */
  pp->ka    = m_ka;
  pp->ra    = m_ra;
  pp->m_max = m_max;
  pp->m_sum = m_sum;
  pp->m_sqr = m_sqr;
  pp->m_len = m_len;
  pp->i_max = i_max;
  pp->i_sum = i_sum;
  pp->i_sqr = i_sqr;
  pp->i_len = i_len;

  pp->m_rs1 = m_rs1;
  pp->pscore = pscore;

  /* record payload in /M/ */
  M->pp = pp;
  M->pp_free = &pp_free;

  return 0;
}
