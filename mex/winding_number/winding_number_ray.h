#ifndef IGL_WINDING_NUMBER_RAY_H
#define IGL_WINDING_NUMBER_RAY_H

#include <Eigen/Core>
// Forward declaration
namespace igl
{
  class EmbreeIntersector;
}

#ifdef WITH_CGAL
#include "CGAL_typedefs.h"
#endif

namespace igl
{
  // Compute winding number exactly|approximatively for a closed|non-closed input
  // mesh using a *single* ray.
  //
  // Inputs:
  //   ei  Embree Intersector, previously initialized to hold mesh (V,F)
  //   N  #F by 3 per-face normals corresponding to (V,F), may be unnormalized
  //   p  3d evaluation point (aka ray origin)
  //   dir  3d ray direction
  // Outputs:
  //   w  winding number, 0: outside, 4*pi: once inside
  //   num_rays_shot  total number of rays shot (for marching)
  //
  // Known Bug: Nearly coplanar facets are sometimes missed. Causing this to be
  // off by a factor of ±#misses*1*4*pi. Also if hits occur near edges/vertices
  // then this may cause false positives or misses.
  //
  void winding_number_ray(
    const igl::EmbreeIntersector & ei,
    const Eigen::MatrixXd & N,
    const Eigen::Vector3d & p, 
    const Eigen::Vector3d & dir,
    double & w,
    int & num_rays_shot);
  // Similar but instead of counting the total number of signed hits, count the
  // number of unsigned hits mod 2
  void winding_number_ray_parity(
    const igl::EmbreeIntersector & ei,
    const Eigen::MatrixXd & N,
    const Eigen::Vector3d & p, 
    const Eigen::Vector3d & dir,
    double & w,
    int & num_rays_shot);
  
  #ifdef WITH_CGAL
  // Same as above but using CGAL
  // Inputs:
  //   tree  AABB Tree containing list of triangles tlist
  //   tlist  list of triangles of (V,F)
  //   p  3d evaluation point (aka ray origin)
  //   dir  3d ray direction
  void winding_number_ray(
    const TriTree & tree,
    const TriangleList & tlist,
    const Eigen::Vector3d & p, 
    const Eigen::Vector3d & dir,
    double & w,
    int & num_rays_shot);
  #endif
  
  
  // Compute the winding number approximatively at a given point using *many*
  // rays.
  //
  // Inputs:
  //   ei  Embree Intersector, previously initialized to hold mesh (V,F)
  //   N  #F by 3 per-face normals corresponding to (V,F), may be unnormalized
  //   p  3d evaluation point (aka ray origin)
  //   dir  #rays by 3 list of ray directions
  // Outputs:
  //   w  winding number, 0: outside, 4*pi: once inside
  //   num_rays_shot  total number of rays shot (for marching), should be ≥#rays
  //
  void winding_number_ray(
    const igl::EmbreeIntersector & ei,
    const Eigen::MatrixXd & N,
    const Eigen::Vector3d & p, 
    const Eigen::MatrixXd & D,
    double & w,
    int & num_rays_shot);
  #ifdef WITH_CGAL
  void winding_number_ray(
    const TriTree & tree,
    const TriangleList & tlist,
    const Eigen::Vector3d & p, 
    const Eigen::MatrixXd & D,
    double & w,
    int & num_rays_shot);
  #endif
}

#endif
