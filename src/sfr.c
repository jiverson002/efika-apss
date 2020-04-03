/* SPDX-License-Identifier: MIT */
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "efika/core.h"
#include "efika/impl.h"

#include "efika/core/gc.h"
#include "efika/core/rename.h"
#include "efika/impl/export.h"
#include "efika/impl/rename.h"

#define BRUTE_FORCE_THRESHOLD 256

typedef struct scratchspace
{
  bool *subp;
} ScratchSpace;

typedef struct hyperplane
{
  val_t l;
  ind_t ln;
  ind_t rn;
  ind_t dim;
} Hyperplane;


static inline val_t
dist(Matrix const * const M, ind_t const p0, ind_t const p1)
{
  /* unpack /M/ */
  ind_t const nc = M->nc;
  ind_t const * const ia = M->ia;
  val_t const * const a  = M->a;

  val_t sum = 0.0;
  for (ind_t i = 0; i < nc; i++) {
    val_t const d = a[ia[p0] + i] - a[ia[p1] + i];
    sum += d * d;
  }
  return sum;
}

static inline
Hyperplane select_hyperplane(Matrix const * const I)
{
  /* unpack /I/ */
  ind_t const nr = I->nr;
  ind_t const * const ia = I->ia;
  val_t const * const a  = I->a;

  ind_t const hp = (ind_t)rand() % nr;
  ind_t const n  = ia[hp + 1] - ia[hp];
  ind_t const cn = 1 + ((n - 1) / 2);
  return (Hyperplane){ a[ia[hp] + cn - 1], cn, n - cn, hp };
}

static int
sfrkd_2(
  val_t const minsim,
  Matrix const * const M,
  Matrix const * const I,
  Hyperplane const HP,
  ScratchSpace const SS,
  Vector * const A
);

static int
sfr1d_1(
  val_t const minsim,
  Matrix const * const M,
  Matrix const * const I,
  Vector * const A
)
{
  /* unpack /I/ */
  ind_t const * const ia = I->ia;
  ind_t const * const ja = I->ja;
  val_t const * const a  = I->a;

  if (0 == ia[1] - ia[0])
    return 0;

  val_t const minsim2 = minsim * minsim;

  for (ind_t j = ia[0]; j < ia[1] - 1; j++) {
    for (ind_t k = j + 1; k < ia[1]; k++) {
      /* check that V_k is no more than minsim away from V_j, if so, then we
       * have visited all possible points within the radius of V_j. */
      val_t const v = a[k] - a[j];
      if (v * v > minsim2)
        break;

      /* compute distance between points V_j and V_k and record as a solution,
       * if the distance is less than minsim */
      val_t const djk = dist(M, ja[j], ja[k]);
      if (djk <= minsim2)
        vector_push_back(A, (Solution){ ja[j], ja[k], djk });
    }
  }

  return 0;
}

static int
sfrkd_1(
  val_t const minsim,
  Matrix const * const M,
  Matrix const * const I,
  ScratchSpace const SS,
  Vector * const A
)
{
  /*==========================================================================*/
  GC_func_init();
  /*==========================================================================*/

  int ret;

  /* unpack /I/ */
  ind_t const nr = I->nr;
  ind_t const nc = I->nc;
  ind_t const * const ia = I->ia;
  ind_t const * const ja = I->ja;
  val_t const * const a  = I->a;

  /* unpack /SS/ */
  bool * const subp = SS.subp;

  /* return the minimal distance found by the brute-force algorithm */
  if (1 == nr || ia[1] - ia[0] < BRUTE_FORCE_THRESHOLD)
    return sfr1d_1(minsim, M, I, A);

  /* select the hyperplane that determines the left and right sub-problems */
  Hyperplane const HP = select_hyperplane(I);
  ind_t const ln  = HP.ln;
  ind_t const rn  = HP.rn;
  ind_t const dim = HP.dim;

  /* FIXME: This memory will be lost if any GC failure happens. */
  Matrix Il, Ir;
  ret = Matrix_init(&Il);
  GC_assert(!ret);
  ret = Matrix_init(&Ir);
  GC_assert(!ret);

  Il.nr = nr;
  Il.nc = nc;
  Il.ia = GC_malloc((nr + 1) * sizeof(*Il.ia));
  Il.ja = GC_malloc(nr * ln * sizeof(*Il.ja));
  Il.a  = GC_malloc(nr * ln * sizeof(*Il.a));

  Ir.nr = nr;
  Ir.nc = nc;
  Ir.ia = GC_malloc((nr + 1) * sizeof(*Ir.ia));
  Ir.ja = GC_malloc(nr * rn * sizeof(*Ir.ja));
  Ir.a  = GC_malloc(nr * rn * sizeof(*Ir.a));

  /* mark points that will be in left and right view */
  for (ind_t i = 0, j = ia[dim]; j < ia[dim + 1]; i++, j++)
    subp[ja[j]] = i >= ln;

  Il.ia[0] = 0;
  Ir.ia[0] = 0;
  for (ind_t i = 0, sln = 0, srn = 0; i < nr; i++) {
    for (ind_t j = ia[i]; j < ia[i + 1]; j++) {
      if (subp[ja[j]]) {
        Ir.ja[srn] = ja[j];
        Ir.a[srn++] = a[j];
      } else {
        Il.ja[sln] = ja[j];
        Il.a[sln++] = a[j];
      }
    }
    Il.ia[i + 1] = sln;
    Ir.ia[i + 1] = srn;
  }

  /* compute all fixed-radius pairs in the two sub-problems */
  ret = sfrkd_1(minsim, M, &Il, SS, A);
  GC_assert(!ret);
  ret = sfrkd_1(minsim, M, &Ir, SS, A);
  GC_assert(!ret);

  /* compute all fixed-radius pairs in the cross-hyperplane slab */
  ret = sfrkd_2(minsim, M, I, HP, SS, A);
  GC_assert(!ret);

  /* free memory */
  Matrix_free(&Il);
  Matrix_free(&Ir);

  return 0;
}

static int
sfrkd_2(
  val_t const minsim,
  Matrix const * const M,
  Matrix const * const I,
  Hyperplane const HP,
  ScratchSpace const SS,
  Vector * const A
)
{
  /*==========================================================================*/
  GC_func_init();
  /*==========================================================================*/

  int ret;

  /* unpack /I/ */
  ind_t const nr = I->nr;
  ind_t const nc = I->nc;
  ind_t const * const ia = I->ia;
  ind_t const * const ja = I->ja;
  val_t const * const a  = I->a;

  /* unpack /HP/ */
  val_t const l   = HP.l ;
  ind_t const ln  = HP.ln;
  ind_t const dim = HP.dim;

  /* unpack /SS/ */
  bool * const subp = SS.subp;

  /* mark all the points for which |x-l| <= t */
  ind_t sn = 0;
  for (ind_t j = ia[dim]; j < ia[dim + 1]; j++) {
    subp[ja[j]] = fabsf(a[j] - l) <= minsim;
    sn += subp[ja[j]];
  }

  /* return now if there will be no points in the slab */
  if (!sn)
    return 0;

  /* allocate memory */
  ind_t * const T = GC_malloc((nr - 1) * sizeof(*T));

  Matrix II;
  ret = Matrix_init(&II);
  GC_assert(!ret);

  II.nr = nr - 1;
  II.nc = nc;
  II.ia = GC_malloc(nr * sizeof(*II.ia));
  II.ja = GC_malloc((nr - 1) * sn * sizeof(*II.ja));
  II.a  = GC_malloc((nr - 1) * sn * sizeof(*II.a));

  II.ia[0] = 0;
  for (ind_t i = 0, ni = 0, sln = 0; i < nr; i++) {
    if (i == dim)
      continue;

    for (ind_t j = ia[i]; j < ia[i + 1]; j++) {
      if (subp[ja[j]]) {
        II.ja[sln] = ja[j];
        II.a[sln++] = a[j];
      }
    }

    II.ia[++ni] = sln;
  }

  /* record current number of solutions, so "false-drops" can be reconciled
   * later. */
  size_t size = A->size;

  /* compute all fixed-radius pairs in the slab sub-problem */
  ret = sfrkd_1(minsim, M, &II, SS, A);
  GC_assert(!ret);

  /* re-mark the left and right sub-problems */
  for (ind_t i = 0, j = ia[dim]; j < ia[dim + 1]; i++, j++)
    subp[ja[j]] = i >= ln;

  /* remove any "false-drops" */
  for (size_t i = size; i < A->size; i++) {
    Solution s = vector_get(A, i);
    if (subp[s.p0] != subp[s.p1])
      A->mem[size++] = s;
  }
  A->size = size;

  /* free memory */
  Matrix_free(&II);

  return 0;
}


EFIKA_IMPL_EXPORT int
Impl_sfr1d(val_t const minsim, Matrix const * const M, Matrix const * const I,
           Vector * const A)
{
  return sfr1d_1(minsim, M, I, A);
}

EFIKA_IMPL_EXPORT int
Impl_sfrkd(val_t const minsim, Matrix const * const M, Matrix const * const I,
           Vector * const A)
{
  /*==========================================================================*/
  GC_func_init();
  /*==========================================================================*/

  bool * const subp = GC_malloc(M->nr * sizeof(*subp));

  int err = sfrkd_1(minsim, M, I, (ScratchSpace){ subp }, A);
  GC_assert(!err);

  return 0;
}
