#include "helmFE.h"
#include <Eigen/Eigen>
#include <Eigen/CholmodSupport>
#include <iostream>

namespace helm {
// Standard function is:
// x
// y
// 1 - x - y

Eigen::Matrix3<Scalar> standard_pdv{}; // TODO: fill the element.

Eigen::Matrix<Scalar, 2, 3> std_dpdxy{{1, 0, -1}, {0, 1, -1}};

Eigen::Matrix3<Scalar> std_pv{{1.0 / 12.0, 1.0 / 24, 1.0 / 24},
                              {1.0 / 24, 1.0 / 12, 1.0 / 24},
                              {1.0 / 24, 1.0 / 24, 1.0 / 12}};

FiniteElement &FiniteElement::setMesh(const Mesh &mesh) {
  mesh_ = mesh;
  return *this;
}

void FiniteElement::solve() {
  stf_coefs_.clear();
  buildStiffness();
  buildLoad();
  buildRobin();
  buildDirichlet();
  stiffness_.resize(mesh_.position_.size(), mesh_.position_.size());
  stiffness_.setFromTriplets(stf_coefs_.begin(), stf_coefs_.end());
  std::cerr << "Solve Eq." << std::endl;
  Eigen::CholmodSimplicialLDLT<SparseMatrix> solver;
  solver.compute(stiffness_);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Fail: " << solver.info() << std::endl;
  }
  values_ = solver.solve(load_);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solve Fail: " << solver.info() << std::endl;
  }
}

void FiniteElement::buildStiffness() {
  for (auto t : mesh_.triangles_) {
    auto pi = mesh_.position_[t(0)];
    auto pj = mesh_.position_[t(1)];
    auto pk = mesh_.position_[t(2)];
    Eigen::Matrix3<Scalar> pijk = Eigen::Matrix3<Scalar>::Ones();
    pijk.block<2, 1>(1, 0) = pi;
    pijk.block<2, 1>(1, 1) = pj;
    pijk.block<2, 1>(1, 2) = pk;
    Scalar area2 = (pijk.determinant());

    // First part: << Dp , Dv >> Inner product
    Scalar dxdr = (pj - pk).y() / area2;
    Scalar dydr = (pk - pi).y() / area2;
    Scalar dxdz = (pk - pj).x() / area2;
    Scalar dydz = (pi - pk).x() / area2;
    Eigen::Matrix2<Scalar> dxy_drz{{dxdr, dydr}, {dxdz, dydz}};
    Eigen::Matrix<Scalar, 2, 3> dpi_drz = dxy_drz * std_dpdxy;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        // Integrate[ Dp dot Dv ]
        Scalar ip_local = dpi_drz.col(i).dot(dpi_drz.col(j));
        Scalar val = ip_local * 0.5 * abs(area2);
        // TODO: Not tested.
        // Second part: Integrate [ 1/r dpdr * v ]
        Scalar dpdr = dpi_drz(0, i);
        Scalar int_dpdr_v = abs(area2) / Scalar(6);
        val += -int_dpdr_v * dpdr / (1e-7 + pi.x());

        // Third part: Integrate [k^2 p * v]
        Scalar k2 = equation_->k2(mesh_.position_[t(j)]);
        Scalar int_p_v = std_pv(j, i) * abs(area2);
        val += -k2 * int_p_v;
        Triplet coef{t(j), t(i), val};
        stf_coefs_.push_back(coef);
      }
    }
  }
}

void FiniteElement::buildLoad() {
  load_.resize(mesh_.position_.size());
  load_.setZero();

  for (auto t : mesh_.triangles_) {
    auto pi = mesh_.position_[t(0)];
    auto pj = mesh_.position_[t(1)];
    auto pk = mesh_.position_[t(2)];
    Eigen::Matrix3<Scalar> pijk = Eigen::Matrix3<Scalar>::Ones();
    pijk.block<2, 1>(1, 0) = pi;
    pijk.block<2, 1>(1, 1) = pj;
    pijk.block<2, 1>(1, 2) = pk;
    Scalar area2 = abs(pijk.determinant());
    for (int i = 0; i < 3; ++i) {
      auto pi = mesh_.position_[t(i)];
      load_(t(i)) += area2 * equation_->f(pi) / 6;
    }
  }
}

void FiniteElement::buildDirichlet() {
  for (auto vid : dirichlet_bc_->vert_ids_) {
    Vec2 p = mesh_.position_[vid];
    Scalar v = dirichlet_bc_->f(p);
    for (auto &trip : stf_coefs_) {
      if (trip.col() == vid) {
        load_(trip.row()) -= v * trip.value();
      }

      if (trip.row() == vid || trip.col() == vid) {
        trip = Triplet(trip.row(), trip.col(), 0);
      }
    }
    stf_coefs_.push_back(Triplet(vid, vid, 1));
    load_(vid) = v;
  }
}

void FiniteElement::buildRobin() {
  // TODO:
}

} // namespace helm
