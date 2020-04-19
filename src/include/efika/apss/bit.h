/* SPDX-License-Identifier: MIT */
#ifndef EFIKA_APSS_BIT_H
#define EFIKA_APSS_BIT_H 1

#define BIT_lsb(i) ((i) & -(i))

#define BIT_sum(n, a, i, op, out)\
  do {\
    ind_t i_ = i + 1;\
    val_t res_ = 0.0;\
    while (i_ > 0)\
      res_ = op(a[i_], res_), i_ -= BIT_lsb(i_);\
    out = res_;\
  } while (0)

#define BIT_add(n, a, i, k, op)\
  do {\
    ind_t i_ = i + 1;\
    val_t k_ = k;\
    while (i_ <= n)\
      a[i_] = op(k_, a[i_]), i_ += BIT_lsb(i_);\
  } while (0)

#endif /* EFIKA_APSS_BIT_H */
