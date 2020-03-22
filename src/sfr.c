/* SPDX-License-Identifier: MIT */
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "efika/impl/sfr.h"

#define BRUTE_FORCE_THRESHOLD 256

typedef unsigned const * View;

typedef struct subproblem {
  unsigned k;
  unsigned n;
  unsigned *m;
  View const *V;
} Subproblem;

typedef struct scratchspace {
  bool *subp;
  View * scratch;
} ScratchSpace;

typedef struct hyperplane {
  float l;
  unsigned ln;
  unsigned rn;
  unsigned dim;
} Hyperplane;


struct kv {
  float k;
  unsigned v;
};

static int kvcmp(void const * const ap, void const * const bp) {
  return ((struct kv const*)ap)->k < ((struct kv const*)bp)->k ? -1 : 1;
}


static inline float dist(Problem const P, unsigned const p0, unsigned const p1) {
  unsigned const k            = P.k;
  float const (*const mem)[5] = P.mem;

  float sum = 0.0;
  for (unsigned i = 0; i < k; i++) {
    float const d = mem[p0][i] - mem[p1][i];
    sum += d * d;
  }
  return sum;
}

static inline Hyperplane select_hyperplane(Problem const P, Subproblem const SP) {
  float const (*const mem)[5] = P.mem;
  unsigned const n            = SP.n;
  unsigned const * const m    = SP.m;
  View const * const V        = SP.V;

  /* TODO: Implement efficient hyperplane selection algorithm */
  //unsigned const hp = 0;
  unsigned const hp = (unsigned)rand() % SP.k;
  unsigned const cn = 1 + ((n - 1) / 2);
  return (Hyperplane){ mem[V[hp][cn - 1]][m[hp]], cn, n - cn, hp };
}

static void sfrkd_2(Hyperplane const HP, Problem const P, Subproblem const SP, ScratchSpace const SS, Vector * const A);

static void sfr1d_1(Problem const P, Subproblem const SP, Vector * const A) {
  /* unpack... */
  float const t               = P.t;
  float const (*const mem)[5] = P.mem;
  unsigned const n            = SP.n;
  unsigned const * const m    = SP.m;
  View const * const V        = SP.V;

  if (0 == n)
    return;

  float const t2 = t * t;

  for (unsigned i = 0; i < n - 1; i++) {
    for (unsigned j = i + 1; j < n; j++) {
      /* check that V_j is no more than P.t away from V_i, if so, then we have
       * visited all possible points within the radius of V_i. */
      float const v = mem[V[0][j]][m[0]] - mem[V[0][i]][m[0]];
      if (v * v > t2)
        break;

      /* compute distance between points V_i and V_j and record as a solution,
       * if the distance is less than P.t */
      float const dij = dist(P, V[0][i], V[0][j]);
      if (dij <= t2)
        vector_push_back(A, (Solution){ V[0][i], V[0][j], dij });
    }
  }
}

static void sfrkd_1(Problem const P, Subproblem const SP, ScratchSpace const SS, Vector * const A) {
  /* unpack... */
  unsigned const k     = SP.k;
  unsigned const n     = SP.n;
  View const * const V = SP.V;
  bool * const subp    = SS.subp;

  /* return the minimal distance found by the brute-force algorithm */
  if (1 == k || n < BRUTE_FORCE_THRESHOLD) {
    sfr1d_1(P, SP, A);
    return;
  }

  /* select the hyperplane that determines the left and right sub-problems */
  Hyperplane const HP = select_hyperplane(P, SP);
  unsigned const ln  = HP.ln;
  unsigned const rn  = HP.rn;
  unsigned const dim = HP.dim;

  /* allocate memory for array of new views */
  View * const Vk = malloc(2 * k * sizeof(*Vk));
  assert(Vk);

  /* allocate memory for new views */
  unsigned * const Vn = malloc((k - 1) * n * sizeof(*Vn));
  assert(Vn);

  /* create view arrays Vl and Vr */
  View * const Vl = Vk;
  View * const Vr = Vk + k;

  /* mark points that will be in left and right view */
  for (unsigned i = 0; i < n; i++)
    subp[V[dim][i]] = i < ln;

  for (unsigned i = 0, j = 0; i < k; i++) {
    unsigned *vn;

    if (i == dim) {
      vn = (unsigned*)V[dim];
    } else {
      vn = Vn + n * j++; /* offset into block allocation */

      unsigned *vl = vn;
      unsigned *vr = vn + ln;
      for (unsigned l = 0; l < n; l++) {
        if (subp[V[i][l]])
          *vl++ = V[i][l];
        else
          *vr++ = V[i][l];
      }
    }

    /* create views vl and vr */
    Vl[i] = vn;
    Vr[i] = vn + ln;
  }

  /* compute all fixed-radius pairs in the two sub-problems */
  sfrkd_1(P, (Subproblem){ k, ln, SP.m, Vl }, SS, A);
  sfrkd_1(P, (Subproblem){ k, rn, SP.m, Vr }, SS, A);

  /* compute all fixed-radius pairs in the cross-hyperplane slab */
  sfrkd_2(HP, P, SP, (ScratchSpace){ subp, Vl }, A);

  /* free memory */
  free(Vn);
  free(Vk);
}

static void sfrkd_2(Hyperplane const HP, Problem const P, Subproblem const SP, ScratchSpace const SS, Vector * const A) {
  /* unpack... */
  float const l               = HP.l ;
  unsigned const ln           = HP.ln;
  unsigned const dim          = HP.dim;
  float const t               = P.t;
  float const (*const mem)[5] = P.mem;
  unsigned const k            = SP.k;
  unsigned const n            = SP.n;
  unsigned const * const m    = SP.m;
  View const * const V        = SP.V;
  bool * const subp           = SS.subp;
  View * const scratch        = SS.scratch;

  /* mark all the points for which |x-l| <= t */
  unsigned sn = 0;
  for (unsigned i = 0; i < n; i++) {
    subp[V[dim][i]] = fabsf(mem[V[dim][i]][m[dim]] - l) <= t;
    sn += subp[V[dim][i]];
  }

  /* return now if there will be no points in the slab */
  if (!sn)
    return;

  unsigned * const M = malloc((k - 1) * sizeof(*M));
  assert(M);
  View * const S = malloc((k - 1) * sizeof(*S));
  assert(S);
  for (unsigned i = 0, j = 0; i < k; i++) {
    if (i == dim)
      continue;

    /* create view */
    unsigned *si = (unsigned*)scratch[i];
    for (unsigned ii = 0; ii < n; ii++)
      if (subp[V[i][ii]])
        *si++ = V[i][ii];

    M[j]   = m[i];
    S[j++] = scratch[i];
  }

  /* record current number of solutions, so "false-drops" can be reconciled
   * later. */
  size_t size = A->size;

  /* compute all fixed-radius pairs in the slab sub-problem */
  sfrkd_1(P, (Subproblem){ k - 1, sn, M, S }, SS, A);

  /* re-mark the left and right sub-problems */
  for (unsigned i = 0; i < n; i++)
    subp[V[dim][i]] = i < ln;

  /* remove any "false-drops" */
  for (size_t i = size; i < A->size; i++) {
    Solution s = vector_get(A, i);
    if (subp[s.p0] != subp[s.p1])
      A->mem[size++] = s;
  }
  A->size = size;

  /* free memory */
  free(M);
  free(S);
}


static Subproblem make_subproblem(Problem const P) {
  /* allocate key/value array */
  struct kv * const kv = malloc(P.n * sizeof(*kv));
  assert(kv);

  View * const Vk = malloc(P.k * sizeof(*Vk));
  assert(Vk);

  unsigned * const Vn = malloc(P.k * P.n * sizeof(*Vn));
  assert(Vn);

  unsigned * const m = malloc(P.k * sizeof(*m));
  assert(m);

  for (unsigned i = 0; i < P.k; i++) {
    /* sort points according to this dimension */
    for (unsigned j = 0; j < P.n; j++)
      kv[j].k = P.mem[j][i], kv[j].v = j;
    qsort(kv, P.n, sizeof(*kv), kvcmp);

    /* create view for this dimension */
    unsigned *vn = Vn + i * P.n;

    /* create this dimension */
    for (unsigned j = 0; j < P.n; j++)
      vn[j] = kv[j].v;

    /* record view for this dimension */
    Vk[i] = vn;
    m[i]  = i;
  }

  free(kv);

  return (Subproblem){ P.k, P.n, m, Vk };
}

static ScratchSpace make_scratchspace(Problem const P) {
  bool * const subp = malloc(P.n * sizeof(*subp));
  assert(subp);

  return (ScratchSpace){ subp, NULL };
}

void sfr1d(Problem const P, Vector * const A) {
  Subproblem const SP = make_subproblem(P);

  sfr1d_1(P, SP, A);

  free((void*)SP.V[0]);
  free((void*)SP.V);
}

void sfrkd(Problem const P, Vector * const A) {
  Subproblem const SP = make_subproblem(P);
  ScratchSpace const SS = make_scratchspace(P);

  sfrkd_1(P, SP, SS, A);

  free(SP.m);
  free((void*)SP.V[0]);
  free((void*)SP.V);
  free(SS.subp);
}
