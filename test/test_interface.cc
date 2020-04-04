// SPDX-License-Identifier: MIT
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <new>
#include <random>
#include <stdexcept>

#include <gtest/gtest.h>

#include "efika/data.h"
#include "efika/core.h"
#include "efika/impl.h"
#include "efika/io.h"

namespace impl {

struct sfr1d {
  int operator()(EFIKA_val_t const minsim, EFIKA_Matrix const * const M,
                 EFIKA_Matrix const * const I, Vector * const A) {
    return EFIKA_Impl_sfr1d(minsim, M, I, A);
  }
};

struct sfrkd {
  int operator()(EFIKA_val_t const minsim, EFIKA_Matrix const * const M,
                 EFIKA_Matrix const * const I, Vector * const A) {
    return EFIKA_Impl_sfrkd(minsim, M, I, A);
  }
};

} // impl

namespace {

template <typename TypeParam>
class interface : public ::testing::Test {
  public:
    interface(float const t, const std::string &f)
      : minsim_(t), filename_(f) { }

    void SetUp() override {
      int err;

      err = EFIKA_Matrix_init(&M_);
      if (err)
        throw std::runtime_error("Could not initialize matrix");
      err = EFIKA_Matrix_init(&I_);
      if (err)
        throw std::runtime_error("Could not initialize matrix");

      FILE * fp = fopen((filename_).c_str(), "r");
      if (!fp)
        throw std::invalid_argument("Cannot open `" + filename_ + "' for reading");

      err = EFIKA_IO_cluto_load(fp, &M_);
      if (err)
        throw std::runtime_error("Could not load `" + filename_ + "'");

      fclose(fp);

      err = EFIKA_Matrix_sort(&M_, EFIKA_ASC | EFIKA_COL);
      if (err)
        throw std::runtime_error("Could not sort matrix");

      err = EFIKA_Matrix_iidx(&M_, &I_);
      if (err)
        throw std::runtime_error("Could not transpose matrix");

      err = EFIKA_Matrix_sort(&I_, EFIKA_ASC | EFIKA_VAL);
      if (err)
        throw std::runtime_error("Could not sort matrix");
    }

    void TearDown() override {
      EFIKA_Matrix_free(&M_);
      EFIKA_Matrix_free(&I_);
    }

    void TestBody() override {
      int err;

      // declare solution vector
      Vector A;

      // find all fixed-radius pairs using /brute-force/ algorithm
      A = vector_new();
      err = impl::sfr1d()(minsim_, &M_, &I_, &A);
      ASSERT_EQ(err, 0);
      auto const size_brute_force = A.size;
      vector_delete(&A);

      // find all fixed-radius pairs using /efficient/ algorithm
      A = vector_new();
      err = TypeParam()(minsim_, &M_, &I_, &A);
      ASSERT_EQ(err, 0);
      auto const size_efficient = A.size;
      vector_delete(&A);

      ASSERT_EQ(size_brute_force, size_efficient);
    }

  private:
    float minsim_;
    std::string filename_;
    EFIKA_Matrix M_;
    EFIKA_Matrix I_;
};

} // namespace

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  auto ndatasets = sizeof(EFIKA_datasets) / sizeof(*EFIKA_datasets);

  for(decltype(ndatasets) i = 0; i < ndatasets; i++) {
    ::testing::RegisterTest(
      "interface", ("sfrkd." + std::string(EFIKA_datasets[i])).c_str(),
      nullptr, nullptr, __FILE__, __LINE__,
      [=]() -> ::testing::Test* {
        return new interface<impl::sfrkd>(0.10,
            std::string(EFIKA_DATA_PATH) + "/" + std::string(EFIKA_datasets[i]));
      });
  }

  return RUN_ALL_TESTS();
}
