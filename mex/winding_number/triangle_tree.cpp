#include "triangle_tree.h"
#include <list>

void triangle_tree(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  TriTree & tree,
  TriangleList & tlist)
{
  assert(F.cols() == 3);
  tlist.clear();

  // Loop over facets
  for(int f = 0;f<F.rows();f++)
  {
    Point3 a(V(F(f,0),0), V(F(f,0),1), V(F(f,0),2));
    Point3 b(V(F(f,1),0), V(F(f,1),1), V(F(f,1),2));
    Point3 c(V(F(f,2),0), V(F(f,2),1), V(F(f,2),2));
    tlist.push_back(Triangle3( a,b,c));
  }
  // constructs AABB tree
  tree.clear();
  tree.insert(tlist.begin(),tlist.end());
}
