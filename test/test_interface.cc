// SPDX-License-Identifier: MIT
#include <iostream>
#include <stdexcept>

#include <gtest/gtest.h>

#include "efika/apss.h"
#include "efika/core.h"
#include "efika/data.h"
#include "efika/io.h"

namespace {

class interface : public ::testing::Test {
  private:
    using impl_func = int (*)(EFIKA_val_t, EFIKA_Matrix*, EFIKA_Matrix*);

  public:
    interface(float const t, const std::string &f, impl_func impl)
      : minsim_(t), filename_(f), impl_(impl) { }

    void SetUp() override {
      int err;

      err = EFIKA_Matrix_init(&M_);
      if (err)
        throw std::runtime_error("Could not initialize matrix");

      err = EFIKA_IO_cluto_load(filename_.c_str(), &M_);
      if (err)
        throw std::runtime_error("Could not load `" + filename_ + "'");

      err = EFIKA_Matrix_comp(&M_);
      if (err)
        throw std::runtime_error("Could not compact column space");

      err = EFIKA_Matrix_norm(&M_);
      if (err)
        throw std::runtime_error("Could not normalize matrix");
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
      err = EFIKA_apss_idxjoin(minsim_, &M_, &S);
      ASSERT_EQ(err, 0);
      auto const size_baseline = S.nnz;
      EFIKA_Matrix_free(&S);

      // find all fixed-radius pairs using /efficient/ algorithm
      err = impl_(minsim_, &M_, &S);
      ASSERT_EQ(err, 0);
      auto const size_efficient = S.nnz;
      EFIKA_Matrix_free(&S);

      std::cout << size_baseline << " " << size_efficient << std::endl;

      ASSERT_EQ(size_baseline, size_efficient);
    }

  private:
    float minsim_;
    std::string filename_;
    EFIKA_Matrix M_;
    impl_func impl_;
};

} // namespace

#define REGISTER_TEST_SET(impl)\
  do {\
    auto ndatasets = sizeof(EFIKA_datasets) / sizeof(*EFIKA_datasets);\
    for(decltype(ndatasets) i = 0; i < ndatasets; i++) {\
      ::testing::RegisterTest(#impl, EFIKA_datasets[i],\
        nullptr, nullptr, __FILE__, __LINE__,\
        [=]() -> ::testing::Test* {\
          return new interface(0.90,\
            std::string(EFIKA_DATA_PATH"/") + std::string(EFIKA_datasets[i]),\
            EFIKA_apss_ ## impl);\
        });\
    }\
  } while (0)

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  //REGISTER_TEST_SET(allpairs);
  REGISTER_TEST_SET(l2ap);
  REGISTER_TEST_SET(mmjoin);
  REGISTER_TEST_SET(nova);

  return RUN_ALL_TESTS();
}
