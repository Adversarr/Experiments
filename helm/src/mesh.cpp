#include "mesh.h"
namespace helm {
Mesh GenerateUniform(int div_x, int div_y, Scalar x, Scalar y) {
  Mesh mesh;
  mesh.position_.resize((div_x + 1) * (div_y + 1));
  mesh.triangles_.resize(div_x * div_y * 2);
  Scalar dx = x / static_cast<Scalar>(div_x);
  Scalar dy = y / static_cast<Scalar>(div_y);
  for (int i = 0; i <= div_x; ++i) {
    for (int j = 0; j <= div_y; ++j) {
      Scalar xi = static_cast<Scalar>(i) * dx;
      Scalar yi = static_cast<Scalar>(j) * dy;
      mesh.position_[i * (div_y + 1) + j] = Vec2{xi, yi};
    }
  }

  size_t cnt = 0;
  for (int i = 0; i < div_x; ++i) {
    for (int j = 0; j < div_y; ++j) {
      Eigen::Index lt = i * (div_y + 1) + j;
      auto rt = lt + 1;
      auto lb = lt + (div_y + 1);
      auto rb = lb + 1;
      mesh.triangles_[cnt++] = IVec3{lt, lb, rb};
      mesh.triangles_[cnt++] = IVec3{rb, lt, rt};
    }
  }
  return mesh;
}

} // namespace helm
