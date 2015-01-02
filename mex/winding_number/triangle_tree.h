#ifndef TRIANGLE_TREE_H
#define TRIANGLE_TREE_H

#include <Eigen/Core>
#include "CGAL_typedefs.h"

// Construct an AABB tree using CGAL and fill it with the triangles of a
// triangle mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of vertex positions
//   F  #F by 3 list of triangle indices
// Outputs:
//   Tree  AABB tree containing triangles
//   tlist  list of triangles
void triangle_tree(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  TriTree & tree,
  TriangleList & tlist);
#endif
