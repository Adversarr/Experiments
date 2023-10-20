#pragma once
#include "eqn.h"

namespace helm {
class FiniteDifference : public SolverBase {

public:
  using SolverBase::SolverBase;

  void solve() override;

  FiniteDifference &setRDivision(int r) {
    r_div_ = r;
    return *this;
  }

  FiniteDifference &setZDivision(int z) {
    z_div_ = z;
    return *this;
  }

private:
  Vector grid_value_;
  SparseMatrix spm_;

  int r_div_, z_div_;
};
} // namespace helm
