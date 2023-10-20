#include "helmFE.h"
#include <Eigen/Eigen>
#include <iostream>

using namespace helm;
using namespace std;
int N = 200;
Scalar R = 50, Z = 50;

class MyEqn : public HelmEqn {
  virtual Scalar f(Vec2 p) override final {
    // return 0;
    if (abs(p.x()) < 2 * R / N && abs(p.y() - 0.5 * Z) < 2 * Z / N) {
      return -1;
    } else {
      return 0;
    }
  }
};

class MyDir : public DirichletBC {
  virtual Scalar f(Vec2 p) override final {
    if (abs(p.x()) < 2 * R / N && abs(p.y() - 0.5 * Z) < 2 * Z / N) {
      return -1;
    } else {
      return 0;
    }
  }
};

int main(int argc, char **argv) {
  Eigen::initParallel();
  Eigen::setNbThreads(10);
  auto mesh = helm::GenerateUniform(N, N, R, Z);
  helm::FiniteElement fe;
  MyEqn eqn;
  MyDir dir;
  helm::RobinBC robin;
  // Setup the dirichlet boundaries.
  for (auto id = 0; id < mesh.position_.size(); ++id) {
    auto p = mesh.position_[id];
    if ( // (abs(p.x()) < 2 * R / N && abs(p.y() - 0.5 * Z) < 2 * Z / N) ||
        p.y() < Z / N  // || p.y() >= Z - Z / N
    ) {
      dir.vert_ids_.push_back(id);
    }
  }
  fe.setMesh(mesh);
  fe.equation_ = &eqn;
  fe.dirichlet_bc_ = &dir;
  fe.robin_bc_ = &robin;
  fe.solve();

  Eigen::Matrix<helm::Scalar, Eigen::Dynamic, Eigen::Dynamic> result =
      fe.values_;

  result.resize(N + 1, N + 1);
  // std::cout << fe.stiffness_ << std::endl;
  std::cout << result << std::endl;
  return 0;
}
