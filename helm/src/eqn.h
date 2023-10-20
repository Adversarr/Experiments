#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

namespace helm {

using Scalar = double;
using Vec2 = Eigen::Vector<Scalar, 2>;

using IVec3 = Eigen::Vector<int, 3>;

// 2 types of BC are supported: Dirichlet and Robin
// (Neumann is included in Robin bc. see below)
inline Scalar sqr(Scalar x) noexcept { return x * x; }

/**
 * @class DirichletBC
 * @brief Dirichlet boundary condition:
 *          u(x, y) = f
 */
struct DirichletBC {
  /**
   * @brief Returns the value of boundary at given point.
   *
   * @param X
   * @return
   */
  virtual Scalar f(Vec2 X) { return 0; }
  std::vector<Eigen::Index> vert_ids_;
};

/**
 * @class RobinBC
 * @brief Robin bc:
 *      diff u
 *      ------ + b(x, y) u = g(x, y)
 *      diff n
 */
struct RobinBC {
  /**
   * @brief Returns the value of {b, g} at given point.
   *
   * @param X
   * @return
   */
  virtual Vec2 f(Vec2 X) { return Vec2::Zero(); }
  std::vector<Eigen::Index> vert_ids_;
};

using HomoDirichletBC = DirichletBC;

using HomoNeumannBC = RobinBC;

/**
 * @class HelmEqn
 * @brief Helmholtz equation:
 *  1   diff u   diff^2 u   diff^2 z
 *  - · ------ + -------- + -------- + k^2 u = f(r, z)
 *  r   diff r   diff r^2   diff z^2
 */
struct HelmEqn {
  virtual Scalar k2(Vec2 X) { return sqr(2 * 3.14159 * 300 / 1500); }
  virtual Scalar f(Vec2 X) {
    auto x = X.x(), y = X.y();
    return -((2 * cos(x) - (-1 + x) * sin(x)) * (y - 1) * sin(y) +
             (x - 1) * sin(x) * (2 * cos(y) - (-1 + y) * sin(y)));
  }
};

using TriangleList = std::vector<IVec3>;
using VertexPosition = std::vector<Vec2>;
using Vector = Eigen::Vector<Scalar, Eigen::Dynamic>;
using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;

class SolverBase {
public:
  virtual void solve() = 0;

  // protected:
  DirichletBC *dirichlet_bc_;
  RobinBC *robin_bc_;
  HelmEqn *equation_;
};
} // namespace helm
