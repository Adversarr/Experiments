#pragma once
#include <Eigen/Core>

template <int dim>
Eigen::Index sub2ind(Eigen::Vector<Eigen::Index, dim> size,
                     Eigen::Vector<Eigen::Index, dim> sub) {
  Eigen::Index ind = sub(0);
  for (Eigen::Index i = 1; i < size.size(); ++i) {
    ind = ind * size(i) + sub(i);
  }
  return ind;
}


