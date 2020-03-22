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

#include "efika/data/path.h"
#include "efika/impl/sfr.h"

namespace impl {

struct sfr1d {
  void operator()(Problem const P, Vector * const A) {
    ::sfr1d(P, A);
  }
};

struct sfrkd {
  void operator()(Problem const P, Vector * const A) {
    ::sfrkd(P, A);
  }
};

} // impl

namespace {

template <typename TypeParam>
class interface : public ::testing::Test {
  public:
    interface(float const t, const std::string &f)
      : threshold_(t), filename_(f) { }

    void SetUp() override {
      unsigned n, k;

      std::ifstream file(filename_);
      if (!file.is_open())
        throw std::invalid_argument("Cannot open `" + filename_ + "' for reading");

      file >> n >> k;
      if (file.fail())
        throw std::invalid_argument("Cannot read `" + filename_);

      mem_ = std::make_unique<float[][5]>(n);

      for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < k; j++) {
          file >> mem_.get()[i][j];
          if (file.fail())
            throw std::invalid_argument("Cannot read `" + filename_);
        }
      }

      P_ = { threshold_, k, n, mem_.get() };
    }

    void TestBody() override {
      // declare solution vector
      Vector A;

      // find all fixed-radius pairs using /brute-force/ algorithm
      A = vector_new();
      impl::sfr1d()(P_, &A);
      auto const size_brute_force = A.size;
      vector_delete(&A);

      // find all fixed-radius pairs using /efficient/ algorithm
      A = vector_new();
      TypeParam()(P_, &A);
      auto const size_efficient = A.size;
      vector_delete(&A);

      ASSERT_EQ(size_brute_force, size_efficient);
    }

  private:
    float threshold_;
    std::string filename_;
    Problem P_;
    std::unique_ptr<float[][5]> mem_;
};

} // namespace

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  auto ndatasets = sizeof(datasets) / sizeof(*datasets);

  for(decltype(ndatasets) i = 0; i < ndatasets; i++) {
    ::testing::RegisterTest(
      "interface", ("sfrkd." + std::string(datasets[i])).c_str(),
      nullptr, nullptr, __FILE__, __LINE__,
      [=]() -> ::testing::Test* {
        return new interface<impl::sfrkd>(0.10,
            std::string(EFIKA_DATA_PATH) + "/" + std::string(datasets[i]));
      });
  }

  return RUN_ALL_TESTS();
}
