#pragma once
#include "eqn.h"
#include "mesh.h"
namespace helm {
class FiniteElement : public SolverBase {
public:
  using Triplet = Eigen::Triplet<Scalar>;
  using TripList = std::vector<Triplet>;
  virtual void solve() override;

  FiniteElement& setMesh(const Mesh& mesh);
// private:
  void buildStiffness();
  void buildLoad();
  void buildDirichlet();
  void buildRobin();

// private:
  SparseMatrix stiffness_;

  TripList stf_coefs_;
  Vector load_;
  Vector values_;
  Mesh mesh_;
};
} // namespace helm
