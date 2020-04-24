/* SPDX-License-Identifier: MIT */
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/nova.h"
#include "efika/apss/export.h"
#include "efika/apss/rename.h"
#include "efika/core/blas.h"
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
/* Precompute prefix statistics for a single row. */
/*----------------------------------------------------------------------------*/
static inline Stat
stat_row(
  val_t const minsim,
  ind_t const n,
  ind_t const * const restrict ja,
  val_t const * const restrict a,
  val_t       * const restrict rowmax,
  val_t       * const restrict rowsum,
  val_t       * const restrict rowlen,
  val_t       * const restrict rowrs1,
  val_t       * const restrict colmax
)
{
  bool have_stat = false;
  Stat stat;
  val_t rs1 = 0.0;
  val_t max = 0.0;
  val_t sum = 0.0;
  val_t sqr = 0.0;
  val_t ub  = 0.0;
  val_t pscore;

  for (ind_t j = 0; j < n; j++) {
    rs1 += a[j] * colmax[ja[j]]; /* Bayardo   */
    max  = max(max, a[j]);       /* Awekar    */
    sum += a[j];                 /* ...       */
    sqr += a[j] * a[j];          /* Anastasiu */
                                 /* ...       */
    rowrs1[j] = rs1;
    rowmax[j] = max;
    rowsum[j] = sum;
    rowlen[j] = sqrtv(sqr);

    colmax[ja[j]] = max(colmax[ja[j]], a[j]);

    /* record pscore before updating ub */
    pscore = ub;

    /* dot-product upper-bound */
    ub = min(
      rowrs1[j], /* Bayardo  */
      rowlen[j]  /* Anastasiu */
    );

    if (!have_stat && ub >= minsim) {
      stat = (Stat){ j, pscore };
      have_stat = true;
    }
  }

  /* return prefix split */
  return have_stat ? stat : (Stat){ n, pscore };
}

/*----------------------------------------------------------------------------*/
/* Precompute prefix for a single row. */
/*----------------------------------------------------------------------------*/
static inline ind_t
split_row(
  val_t const minsim,
  ind_t const n,
  ind_t const * const restrict ja,
  val_t const * const restrict a,
  val_t       * const restrict rowlen,
  val_t       * const restrict colmax
)
{
  val_t sqr = 0.0;
  val_t rs1 = 0.0;

  for (ind_t j = 0; j < n; j++) {
    sqr += a[j] * a[j];
    rs1 += a[j] * colmax[ja[j]];

    rowlen[j] = sqrtv(sqr);

    if ( rowlen[j] >= minsim   /* Bayardo   */
      && rs1       >= minsim ) /* Anastasiu */
    {
      return j;
    }

    colmax[ja[j]] = max(colmax[ja[j]], a[j]);
  }

  /* return prefix split */
  return n - 1;
}

/*----------------------------------------------------------------------------*/
/*! Precompute row prefixes. */
/*----------------------------------------------------------------------------*/
static int
filter(
  val_t const minsim,
  ind_t const nr,
  ind_t const nc,
  ind_t const * const restrict ia,
  ind_t const * const restrict ja,
  val_t const * const restrict a,
  ind_t       * const restrict ka,
  val_t       * const restrict max,
  val_t       * const restrict sum,
  val_t       * const restrict len,
  val_t       * const restrict rs1,
  ind_t       * const restrict psplit,
  val_t       * const restrict pscore
)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* allocate scratch memory */
  val_t * const colmax1 = GC_calloc(nc, sizeof(*colmax1));
  val_t * const colmax2 = GC_calloc(nc, sizeof(*colmax2));

  for (ind_t ip1 = nr; ip1 > 0; ip1--) {
    ind_t i = ip1 - 1;
    ind_t j = nr - ip1;

    /* find prefix split and pscore for each row_i */
    Stat const stat = stat_row(minsim, ia[i + 1] - ia[i], ja + ia[i], a + ia[i],
                               max + ia[i], sum + ia[i], len + ia[i],
                               rs1 + ia[i], colmax1);
    ka[i]     = ia[i] + stat.n;
    pscore[i] = stat.pscore;

    psplit[j] = ia[j] + split_row(minsim, ia[j + 1] - ia[j], ja + ia[j],
                                  a + ia[j], len + ia[j], colmax2);
  }
  /* free scratch memory */
  GC_free(colmax1);
  GC_free(colmax2);

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
  val_t const * const restrict len,
  val_t const * const restrict rs1,
  ind_t       * const restrict ia1,
  ind_t       * const restrict ja1,
  val_t       * const restrict acsc,
  val_t       * const restrict max1,
  val_t       * const restrict sum1,
  val_t       * const restrict len1,
  val_t       * const restrict rs11
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
      len1[ia1[ja[j]]]   = len[j];
      rs11[ia1[ja[j]]++] = rs1[j];
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
  val_t  const * const m_len,
  val_t  const * const m_rs1,
  Matrix       * const I,
  val_t ** const i_max,
  val_t ** const i_sum,
  val_t ** const i_len,
  val_t ** const i_rs1
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
               *i_len = GC_malloc(i_nnz * sizeof(*i_len));
               *i_rs1 = GC_malloc(i_nnz * sizeof(*i_rs1));

  csrcsc(m_nr, m_nc, m_ia, m_ja, m_ka, m_a, m_max, m_sum, m_len, m_rs1, i_ia,
         i_ja, i_a, *i_max, *i_sum, *i_len, *i_rs1);

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
  free(pp->m_len);
  free(pp->m_rs1);

  free(pp->i_max);
  free(pp->i_sum);
  free(pp->i_len);
  free(pp->i_rs1);

  free(pp->psplit);
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
  val_t * const m_max = GC_malloc(m_nnz * sizeof(*m_max));
  val_t * const m_sum = GC_malloc(m_nnz * sizeof(*m_sum));
  val_t * const m_len = GC_malloc(m_nnz * sizeof(*m_len));
  val_t * const m_rs1 = GC_malloc(m_nnz * sizeof(*m_rs1));
  ind_t * const psplit = GC_malloc(m_nr * sizeof(*psplit));
  val_t * const pscore = GC_malloc(m_nr * sizeof(*pscore));
  val_t * i_max;
  val_t * i_sum;
  val_t * i_len;
  val_t * i_rs1;

  /* init /I/ */
  err = Matrix_init(I);
  GC_assert(!err);

  /* precompute row prefixes and their statistics */
  err = filter(minsim, m_nr, m_nc, m_ia, m_ja, m_a, m_ka, m_max, m_sum, m_len,
               m_rs1, psplit, pscore);
  GC_assert(!err);

  /* precompute dynamic inverted index */
  err = fiidx(M, m_ka, m_max, m_sum, m_len, m_rs1, I, &i_max, &i_sum, &i_len,
              &i_rs1);
  GC_assert(!err);

  /* precompute dynamic inverted index */
  err = diidx(minsim, M, m_ra, m_max, I);
  GC_assert(!err);

  /* record info in /pp/ */
  pp->ka    = m_ka;
  pp->ra    = m_ra;
  pp->m_max = m_max;
  pp->m_sum = m_sum;
  pp->m_len = m_len;
  pp->m_rs1 = m_rs1;

  pp->i_max = i_max;
  pp->i_sum = i_sum;
  pp->i_len = i_len;
  pp->i_rs1 = i_rs1;

  pp->psplit = psplit;
  pp->pscore = pscore;

  /* record payload in /M/ */
  M->pp = pp;
  M->pp_free = &pp_free;

  return 0;
}
