#pragma once


#include "eqn.h"
namespace helm {

struct Mesh {
  VertexPosition position_;
  TriangleList triangles_;
};


/**
 * 0, 0 --- x, 0
 *  |        |
 *  |        |
 * y, 0 --- x, y
 */
Mesh GenerateUniform(int div_x, int div_y, Scalar x, Scalar y);

}
