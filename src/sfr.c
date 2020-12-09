/* SPDX-License-Identifier: MIT */
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "efika/apss.h"
#include "efika/core.h"

#include "efika/apss/rename.h"
#include "efika/core/export.h"
#include "efika/core/gc.h"
#include "efika/core/pp.h"
#include "efika/core/rename.h"


#define BRUTE_FORCE_THRESHOLD 256


typedef struct scratchspace
{
  ind_t n;
  bool *subp;
  val_t *spa;
  ind_t *lbl;
  ind_t *zer;
} ScratchSpace;

typedef struct hyperplane
{
  val_t l;
  ind_t ln;
  ind_t rn;
  ind_t dim;
} Hyperplane;


static int
sfrkd_2(val_t const minsim, Matrix const * const M, Matrix const * const I,
        Hyperplane const HP, ScratchSpace const SS, Vector * const A);

static inline val_t
dist(Matrix const * const M, ind_t const i, ind_t const j)
{
  /* unpack /M/ */
  ind_t const * const ia = M->ia;
  ind_t const * const ja = M->ja;
  val_t const * const a  = M->a;

  ind_t ii, jj;
  val_t res = 0.0;

  for (ii = ia[i], jj = ia[j]; ii < ia[i + 1] && jj < ia[j + 1];) {
    val_t ai, aj;

    if (ja[ii] < ja[jj]) {
      ai = a[ii++];
      aj = 0.0;
    } else if (ja[ii] > ja[jj]) {
      ai = 0.0;
      aj = a[jj++];
    } else {
      ai = a[ii++];
      aj = a[jj++];
    }

    val_t const d = ai - aj;
    res += d * d;
  }

  for (; ii < ia[i + 1]; ii++)
    res += a[ii] * a[ii];

  for (; jj < ia[j + 1]; jj++)
    res += a[jj] * a[jj];

  return res;
}

static inline Hyperplane
select_hyperplane(Matrix const * const I, ScratchSpace const SS)
{
  /* unpack /I/ */
  ind_t const nr = I->nr;
  ind_t const * const ia = I->ia;
  val_t const * const a  = I->a;

  /* unpack /SS/ */
  ind_t const n = SS.n;

  ind_t const hp  = (ind_t)rand() % nr;
  ind_t const nnz = ia[hp + 1] - ia[hp];
  ind_t const ln  = 1 + ((n - 1) / 2);
  val_t const l   = ln <= (n - nnz) ? 0.0 : a[ia[hp] + (ln - (n - nnz)) - 1];
  return (Hyperplane){ l, ln, n - ln, hp };
}


static int
sfr0d_1(val_t const minsim, Matrix const * const M, Vector * const A)
{
  /* unpack /M/ */
  ind_t const nr = M->nr;

  /* ... */
  val_t const minsim2 = minsim * minsim;

  if (0 == nr)
    return 0;

  for (ind_t i = 0; i < nr - 1; i++) {
    for (ind_t j = i + 1; j < nr; j++) {
      /* compute distance between points V_j and V_k and record as a solution,
       * if the distance is less than minsim */
      val_t const dij = dist(M, i, j);
      if (dij <= minsim2)
        vector_push_back(A, (Solution){ i, j, dij });
    }
  }

  return 0;
}

#if 1
static int
sfr1d_1(val_t const minsim, Matrix const * const M, Matrix const * const I,
        ScratchSpace const SS, Vector * const A)
{
  /* unpack /SS/ */
  ind_t const n = SS.n;
  ind_t const * const lbl  = SS.lbl;

  /* ... */
  val_t const minsim2 = minsim * minsim;

  if (0 == n)
    return 0;

  for (ind_t j = 0; j < n - 1; j++) {
    ind_t const vj = lbl[j];

    /* compare this zero against other zeros */
    for (ind_t k = j + 1; k < n; k++) {
      ind_t const vk = lbl[k];

      /* compute distance between points V_j and V_k and record as a solution,
       * if the distance is less than minsim */
      val_t const djk = dist(M, vj, vk);
      if (djk <= minsim2)
        vector_push_back(A, (Solution){ vj, vk, djk });
    }
  }

  return 0;

  (void)I;
}
#else
static int
sfr1d_1(val_t const minsim, Matrix const * const M, Matrix const * const I,
        ScratchSpace const SS, Vector * const A)
{
  /* unpack /I/ */
  ind_t const * const i_ia = I->ia;
  ind_t const * const i_ja = I->ja;
  val_t const * const i_a  = I->a;

  /* unpack /SS/ */
  ind_t const n = SS.n;
  bool        * const subp = SS.subp;
  ind_t const * const lbl  = SS.lbl;
  ind_t       * const zer  = SS.zer;

  /* ... */
  val_t const minsim2 = minsim * minsim;

  if (0 == n)
    return 0;

  /* XXX: Zero values will need to be collected and handled separately from the
   * non-zero values. In this case, along with each other, the set of values
   * that they need to be compared against will be the same for them all, so it
   * can be precomputed. */
  for (ind_t j = i_ia[0]; j < i_ia[1]; j++)
    subp[i_ja[j]] = true;

  ind_t zn = 0;
  for (ind_t i = 0; i < n; i++)
    if (!subp[lbl[i]])
      zer[zn++] = lbl[i];

  for (ind_t j = i_ia[0]; j < i_ia[1]; j++)
    subp[i_ja[j]] = false;

  if (0 == zn)
    goto NOZEROS;

  ind_t kend;
  for (kend = i_ia[0]; kend < i_ia[1]; kend++) {
    val_t const v = i_a[kend];
    if (v * v > minsim2)
      break;
  }

  for (ind_t j = 0; j < zn - 1; j++) {
    ind_t const vj = zer[j];

    /* compare this zero against other zeros */
    for (ind_t k = j + 1; k < zn; k++) {
      ind_t const vk = zer[k];

      /* compute distance between points V_j and V_k and record as a solution,
       * if the distance is less than minsim */
      val_t const djk = dist(M, vj, vk);
      if (djk <= minsim2)
        vector_push_back(A, (Solution){ vj, vk, djk });
    }

    /* compare this zero against non-zeros */
    for (ind_t k = i_ia[0]; k < kend; k++) {
      ind_t const vk = i_ja[k];

      /* compute distance between points V_j and V_k and record as a solution,
       * if the distance is less than minsim */
      val_t const djk = dist(M, vj, vk);
      if (djk <= minsim2)
        vector_push_back(A, (Solution){ vj, vk, djk });
    }
  }

  NOZEROS:
  if (0 == i_ia[1] - i_ia[0])
    return 0;

  for (ind_t j = i_ia[0]; j < i_ia[1] - 1; j++) {
    ind_t const vj = i_ja[j];
    val_t const aj = i_a[j];

    for (ind_t k = j + 1; k < i_ia[1]; k++) {
      /* check that V_k is no more than minsim away from V_j, if so, then we
       * have visited all possible points within the radius of V_j. */
      val_t const v = i_a[k] - aj;
      if (v * v > minsim2)
        break;

      ind_t const vk = i_ja[k];

      /* compute distance between points V_j and V_k and record as a solution,
       * if the distance is less than minsim */
      val_t const djk = dist(M, vj, vk);
      if (djk <= minsim2)
        vector_push_back(A, (Solution){ vj, vk, djk });
    }
  }

  return 0;
}
#endif


static void
mark_1(Hyperplane const HP, Matrix const * const I, ScratchSpace const SS)
{
  /* unpack /HP/ */
  ind_t const ln  = HP.ln;
  ind_t const dim = HP.dim;

  /* unpack /I/ */
  ind_t const * const ia = I->ia;
  ind_t const * const ja = I->ja;

  /* unpack /SS/ */
  ind_t const n = SS.n;
  bool        * const subp = SS.subp;
  ind_t const * const lbl  = SS.lbl;

  /* ... */
  ind_t const nnz = ia[dim + 1] - ia[dim];

  /* mark all non-zeros */
  for (ind_t j = ia[dim]; j < ia[dim + 1]; j++)
    subp[ja[j]] = true;

  if (ln <= (n - nnz)) {
    /* ...all zeros in left... :: so put some zeros in the right */
    ind_t const rnz = (n - nnz) - ln;
    for (ind_t i = 0, j = 0; i < n && j < rnz; i++) {
      if (!subp[lbl[i]]) {
        subp[lbl[i]] = true;
        j++;
      }
    }
  } else {
    /* ...some non-zeros in left... :: so put some non-zeros in the left */
    ind_t const lnnz = ln - (n - nnz);
    for (ind_t j = ia[dim]; j < ia[dim] + lnnz; j++)
      subp[ja[j]] = false;
  }
}

static ind_t
mark_2(val_t const minsim, Hyperplane const HP, Matrix const * const I,
       ScratchSpace const SS)
{
  /* unpack /HP/ */
  val_t const l   = HP.l;
  ind_t const dim = HP.dim;

  /* unpack /I/ */
  ind_t const * const ia = I->ia;
  ind_t const * const ja = I->ja;
  val_t const * const a  = I->a;

  /* unpack /SS/ */
  ind_t const n = SS.n;
  bool        * const subp = SS.subp;
  ind_t const * const lbl  = SS.lbl;

  ind_t sn;
  if (l <= minsim) {
    /* mark everything */
    for (ind_t i = 0; i < n; i++)
      subp[lbl[i]] = true;

    /* ... */
    sn = n;

    /* reset non-zeros that will not be included in the slab */
    for (ind_t j = ia[dim]; j < ia[dim + 1]; j++) {
      if (fabsf(a[j] - l) > minsim) {
        subp[ja[j]] = false;
        sn--;
      }
    }
  } else {
    /* ... */
    sn = 0;

    /* mark only non-zeros that will be included */
    for (ind_t j = ia[dim]; j < ia[dim + 1]; j++) {
      if (fabsf(a[j] - l) <= minsim) {
        subp[ja[j]] = true;
        sn++;
      }
    }
  }

  return sn;
}

static int
efika_extract(Hyperplane const HP, Matrix const * const I,
              ScratchSpace const SS, Matrix * const Il, Matrix * const Ir)
{
  /* ...garbage collected function... */
  GC_func_init();

  int err;

  /* unpack /HP/ */
  ind_t const dim = HP.dim;

  /* unpack /I/ */
  ind_t const nr  = I->nr;
  ind_t const nc  = I->nc;
  ind_t const nnz = I->nnz;
  ind_t const * const ia = I->ia;
  ind_t const * const ja = I->ja;
  val_t const * const a  = I->a;

  /* unpack /SS/ */
  bool const * const subp = SS.subp;

  if (Il) {
    err = Matrix_init(Il);
    GC_assert(!err);
    GC_register_free(Matrix_free, Il);
  }
  err = Matrix_init(Ir);
  GC_assert(!err);
  GC_register_free(Matrix_free, Ir);

  /* pre-count nnz */
  ind_t rnnz = 0;
  for (ind_t i = 0; i < nr; i++) {
    if (!Il && i == dim)
      continue;

    for (ind_t j = ia[i]; j < ia[i + 1]; j++)
      rnnz += subp[ja[j]];
  }
  ind_t lnnz = nnz - rnnz;

  /* allocate memory */
  if (Il) {
    Il->nr  = nr - !Il;
    Il->nc  = nc;
    Il->nnz = lnnz;
    Il->ia  = GC_malloc((nr + 1) * sizeof(*Il->ia));
    Il->ja  = GC_malloc(lnnz * sizeof(*Il->ja));
    Il->a   = GC_malloc(lnnz * sizeof(*Il->a));
  }
  Ir->nr  = nr - !Il;
  Ir->nc  = nc;
  Ir->nnz = rnnz;
  Ir->ia  = GC_malloc((nr + 1) * sizeof(*Ir->ia));
  Ir->ja  = GC_malloc(rnnz * sizeof(*Ir->ja));
  Ir->a   = GC_malloc(rnnz * sizeof(*Ir->a));

  if (Il)
    Il->ia[0] = lnnz = 0;
  Ir->ia[0] = rnnz = 0;

  for (ind_t i = 0, ii = 1; i < nr; i++) {
    if (!Il && i == dim)
      continue;

    for (ind_t j = ia[i]; j < ia[i + 1]; j++) {
      if (subp[ja[j]]) {
        Ir->ja[rnnz]  = ja[j];
        Ir->a[rnnz++] = a[j];
      } else if (Il) {
        Il->ja[lnnz]  = ja[j];
        Il->a[lnnz++] = a[j];
      }
    }

    if (Il)
      Il->ia[ii] = lnnz;
    Ir->ia[ii++] = rnnz;
  }

  return 0;
}

static void
efika_reset(ScratchSpace const SS)
{
  /* unpack /SS/ */
  ind_t const n = SS.n;
  bool        * const subp = SS.subp;
  ind_t const * const lbl  = SS.lbl;

  for (ind_t i = 0; i < n; i++)
    subp[lbl[i]] = false;
}

static void
efika_remap(ScratchSpace const SS, ind_t * const llbl, ind_t * const rlbl)
{
  /* unpack /SS/ */
  ind_t const n = SS.n;
  bool        * const subp = SS.subp;
  ind_t const * const lbl  = SS.lbl;

  /* construct new label mappings */
  for (ind_t i = 0, lj = 0, rj = 0; i < n; i++) {
    if (subp[lbl[i]])
      rlbl[rj++] = lbl[i];
    else if (llbl)
      llbl[lj++] = lbl[i];
  }
}

static size_t
efika_remove(size_t size, ScratchSpace const SS, Vector * const A)
{
  /* unpack /SS/ */
  bool const * const subp = SS.subp;

  for (size_t i = size; i < A->size; i++) {
    Solution s = vector_get(A, i);
    if (subp[s.p0] != subp[s.p1])
      A->mem[size++] = s;
  }

  return size;
}


static int
sfrkd_1(val_t const minsim, Matrix const * const M, Matrix const * const I,
        ScratchSpace const SS, Vector * const A)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* return the minimal distance found by the brute-force algorithm */
  if (1 == I->nr || SS.n < BRUTE_FORCE_THRESHOLD)
    return sfr1d_1(minsim, M, I, SS, A);

  /* select the hyperplane that determines the left and right subproblems */
  Hyperplane const HP = select_hyperplane(I, SS);
  ind_t const ln = HP.ln;
  ind_t const rn = HP.rn;

  /* mark the left and right subproblems */
  mark_1(HP, I, SS);

  /* extract matrix subproblems */
  Matrix Il, Ir;
  int err = efika_extract(HP, I, SS, &Il, &Ir);
  GC_assert(!err);
  GC_register_free(Matrix_free, &Il);
  GC_register_free(Matrix_free, &Ir);

  /* allocate memory for new labels */
  ind_t * const llbl = GC_malloc(ln * sizeof(*llbl));
  ind_t * const rlbl = GC_malloc(rn * sizeof(*rlbl));

  /* construct new label mappings */
  efika_remap(SS, llbl, rlbl);

  /* reset points */
  efika_reset(SS);

  /* compute all fixed-radius pairs in the two subproblems */
  err = sfrkd_1(minsim, M, &Il,
                (ScratchSpace){ ln, SS.subp, SS.spa, llbl, SS.zer }, A);
  GC_assert(!err);
  err = sfrkd_1(minsim, M, &Ir,
                (ScratchSpace){ rn, SS.subp, SS.spa, rlbl, SS.zer }, A);
  GC_assert(!err);

  /* free memory */
  GC_free(llbl);
  GC_free(rlbl);
  GC_free(&Il);
  GC_free(&Ir);

  /* record current number of solutions, so "false-drops" can be reconciled
   * later. */
  size_t size = A->size;

  /* compute all fixed-radius pairs in the cross-hyperplane slab */
  err = sfrkd_2(minsim, M, I, HP, SS, A);
  GC_assert(!err);

  /* mark the left and right subproblems */
  mark_1(HP, I, SS);

  /* remove any "false-drops" */
  A->size = efika_remove(size, SS, A);

  /* reset points */
  efika_reset(SS);

  return 0;
}

static int
sfrkd_2(val_t const minsim, Matrix const * const M, Matrix const * const I,
        Hyperplane const HP, ScratchSpace const SS, Vector * const A)
{
  /* ...garbage collected function... */
  GC_func_init();

  /* compute slab subproblem */
  ind_t const sn = mark_2(minsim, HP, I, SS);

  /* return now if there will be no points in the slab */
  if (!sn) {
    /* reset points */
    efika_reset(SS);

    return 0;
  }

  /* extract matrix subproblem */
  Matrix Is;
  int err = efika_extract(HP, I, SS, NULL, &Is);
  GC_assert(!err);
  GC_register_free(Matrix_free, &Is);

  /* allocate memory for new labels */
  ind_t * const slbl = GC_malloc(sn * sizeof(*slbl));

  /* construct new label mappings */
  efika_remap(SS, NULL, slbl);

  /* reset points */
  efika_reset(SS);

  /* compute all fixed-radius pairs in the slab subproblem */
  err = sfrkd_1(minsim, M, &Is,
                (ScratchSpace){ sn, SS.subp, SS.spa, slbl, SS.zer }, A);
  GC_assert(!err);

  /* free memory */
  GC_free(&Is);
  GC_free(slbl);

  return 0;
}


EFIKA_EXPORT int
apss_sfr0d(val_t const minsim, Matrix const * const M, Vector * const A)
{
  if (!pp_all(M, A))
    return -1;

  return sfr0d_1(minsim, M, A);
}

EFIKA_EXPORT int
apss_sfrkd(val_t const minsim, Matrix const * const M, Matrix const * const I,
           Vector * const A)
{
  /* ...garbage collected function... */
  GC_func_init();

  if (!pp_all(M, I, A))
    return -1;

  bool * const subp = GC_calloc(M->nr, sizeof(*subp));
  val_t * const spa = GC_calloc(M->nc, sizeof(*spa));
  ind_t * const lbl = GC_malloc(M->nr * sizeof(*lbl));
  ind_t * const zer = GC_malloc(M->nr * sizeof(*zer));

  for (ind_t i = 0; i < M->nr; i++)
    lbl[i] = i;

  int err = sfrkd_1(minsim, M, I, (ScratchSpace){ M->nr, subp, spa, lbl, zer },
                    A);
  GC_assert(!err);

  GC_free(subp);
  GC_free(spa);
  GC_free(lbl);
  GC_free(zer);

  return 0;
}
