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

#include "efika/apss.h"
#include "efika/core.h"
#include "efika/data.h"
#include "efika/io.h"

namespace apss {

struct idxjoin {
  static int run(EFIKA_val_t const minsim, EFIKA_Matrix * const M,
                 EFIKA_Matrix * const S)
  {
    return EFIKA_apss_idxjoin(minsim, M, S);
  }
};

struct allpairs {
  static int run(EFIKA_val_t const minsim, EFIKA_Matrix * const M,
                 EFIKA_Matrix * const S)
  {
    return EFIKA_apss_allpairs(minsim, M, S);
  }
};

#if 0
struct sfr {
  int operator()(EFIKA_val_t const minsim, EFIKA_Matrix const * const M,
                 EFIKA_Matrix const * const I, EFIKA_Matrix * const S)
  {
    return EFIKA_apss_sfr(minsim, M, I, S);
  }
};
#endif

} // apss

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

      FILE * fp = fopen((filename_).c_str(), "r");
      if (!fp)
        throw std::invalid_argument("Cannot open `" + filename_ + "' for reading");

      err = EFIKA_IO_cluto_load(fp, &M_);
      if (err)
        throw std::runtime_error("Could not load `" + filename_ + "'");

      err = EFIKA_Matrix_comp(&M_);
      if (err)
        throw std::runtime_error("Could not compact column space");

      err = EFIKA_Matrix_norm(&M_);
      if (err)
        throw std::runtime_error("Could not normalize matrix");

      fclose(fp);
    }

    void TearDown() override {
      EFIKA_Matrix_free(&M_);
    }

    void TestBody() override {
      int err;

      EFIKA_Matrix S;

      err = EFIKA_Matrix_init(&S);
      ASSERT_EQ(err, 0);

      // compute baseline number of pairs
      err = apss::idxjoin::run(minsim_, &M_, &S);
      ASSERT_EQ(err, 0);
      auto const size_idxjoin = S.nnz;
      EFIKA_Matrix_free(&S);

      // find all fixed-radius pairs using /efficient/ algorithm
      err = TypeParam::run(minsim_, &M_, &S);
      ASSERT_EQ(err, 0);
      auto const size_efficient = S.nnz;
      EFIKA_Matrix_free(&S);

      std::cout << size_idxjoin << " " << size_efficient << std::endl;

      ASSERT_EQ(size_idxjoin, size_efficient);
    }

  private:
    float minsim_;
    std::string filename_;
    EFIKA_Matrix M_;
};

} // namespace

#define REGISTER_TEST_SET(impl)\
  do {\
    auto ndatasets = sizeof(EFIKA_datasets) / sizeof(*EFIKA_datasets);\
    for(decltype(ndatasets) i = 0; i < ndatasets; i++) {\
      ::testing::RegisterTest(#impl, EFIKA_datasets[i],\
        nullptr, nullptr, __FILE__, __LINE__,\
        [=]() -> ::testing::Test* {\
          return new interface<apss::impl>(0.90,\
              std::string(EFIKA_DATA_PATH) + "/" + std::string(EFIKA_datasets[i]));\
        });\
    }\
  } while (0)

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  REGISTER_TEST_SET(allpairs);
  //REGISTER_TEST_SET(sfrkd);

  return RUN_ALL_TESTS();
}
