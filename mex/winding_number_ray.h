#ifndef IGL_WINDING_NUMBER_RAY_H
#define IGL_WINDING_NUMBER_RAY_H

#include <Eigen/Core>
// Forward declaration
namespace igl
{
  namespace embree
  {
    class EmbreeIntersector;
  }
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
    const igl::embree::EmbreeIntersector & ei,
    const Eigen::MatrixXd & N,
    const Eigen::Vector3d & p, 
    const Eigen::Vector3d & dir,
    double & w,
    int & num_rays_shot);
  // Similar but instead of counting the total number of signed hits, count the
  // number of unsigned hits mod 2
  void winding_number_ray_parity(
    const igl::embree::EmbreeIntersector & ei,
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
    const igl::embree::EmbreeIntersector & ei,
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
//
// Implementation

#include <igl/embree/EmbreeIntersector.h>
#include <igl/Hit.h>

#include <iterator> 
#include <iostream>

void igl::winding_number_ray_parity(
  const igl::embree::EmbreeIntersector & ei,
  const Eigen::MatrixXd & N,
  const Eigen::Vector3d & p, 
  const Eigen::Vector3d & dir,
  double & w,
  int & num_rays_shot)
{
  using namespace std;
  using namespace Eigen;
  double count = 0.0;
  vector<Hit > hits;
  ei.intersectRay(p.cast<float>(),dir.cast<float>(),hits,num_rays_shot);
  const int m = N.rows();
  int prev_id = -1;
  for(vector<Hit>::iterator hit = hits.begin();
      hit != hits.end();
      hit++)
  {
    const int id = hit->id % m;
#ifdef VERBOSE
    cout<<id+1<<" "<<"("<<hit->u<<" "<<hit->v<<" "<<(1-(hit->u+hit->v))<<")"<<" ";
#endif
    // HACK: skip consequetive hits on same triangle
    if(id == prev_id)
    {
      continue;
    }
    // "entering" or "exiting"
    count++;
    count = ((int)count % 2);
    prev_id = id;
  }
#ifdef VERBOSE
  cout<<endl;
#endif
  w = count;
}

void igl::winding_number_ray(
  const igl::embree::EmbreeIntersector & ei,
  const Eigen::MatrixXd & N,
  const Eigen::Vector3d & p, 
  const Eigen::Vector3d & dir,
  double & w,
  int & num_rays_shot)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  double count = 0.0;
  vector<Hit > hits;
  ei.intersectRay(p.cast<float>(),dir.cast<float>(),hits,num_rays_shot);
  const int m = N.rows();
  int prev_id = -1;
  for(vector<Hit>::iterator hit = hits.begin();
      hit != hits.end();
      hit++)
  {
    const int id = hit->id % m;
#ifdef VERBOSE
    cout<<id+1<<" "<<"("<<hit->u<<" "<<hit->v<<" "<<(1-(hit->u+hit->v))<<")"<<" ";
#endif
    // HACK: skip consequetive hits on same triangle
    if(id == prev_id)
    {
      continue;
    }
    // "entering" or "exiting"
    count += (dir.dot(N.row(id)) > 0 ? 1: -1);
    prev_id = id;
  }
#ifdef VERBOSE
  cout<<endl;
#endif
  w = count;
}


#ifdef WITH_CGAL
void winding_number_ray(
  const TriTree & tree,
  const TriangleList & tlist,
  const Eigen::Vector3d & p, 
  const Eigen::Vector3d & dir,
  double & w,
  int & num_rays_shot)
{
  using namespace std;
  Vector3 d(dir(0),dir(1),dir(2));
  Ray3 ray(Point3(p(0),p(1),p(2)),d);
  // computes all intersected primitives with segment query as primitive ids
  list<Primitive_id> primitives;
  tree.all_intersected_primitives(ray, back_inserter(primitives));
  double count = 0;
  for(
    list<Primitive_id>::const_iterator pit = primitives.begin();
    pit != primitives.end();
    pit++)
  {
    const TriangleListConstIterator tit = *pit;
    //std::size_t index_in_vector = std::distance(tlist.begin(),tit); 
    const Triangle3 & t = *tit;
    // dot product
    double s = d * t.supporting_plane().orthogonal_vector();
    if(s>0)
    {
      count += 1;
    }else
    {
      count += -1;
    }
  }
  w = count;

}
#endif

//void winding_number_ray(
//  const SegTree & tree,
//  const SegmentList & slist,
//  const Eigen::Vector2d & p, 
//  const Eigen::Vector2d & dir,
//  double & w,
//  int & num_rays_shot)
//{
//  using namespace std;
//  using namespace igl;
//  Vector3 d(dir(0),dir(1),dir(2));
//  Ray3 ray(Point3(p(0),p(1),p(2)),d);
//  // computes all intersected primitives with segment query as primitive ids
//  list<Primitive_id> primitives;
//  tree.all_intersected_primitives(ray, back_inserter(primitives));
//  double count = 0;
//  for(
//    list<Primitive_id>::const_iterator pit = primitives.begin();
//    pit != primitives.end();
//    pit++)
//  {
//    const TriangleListConstIterator tit = *pit;
//    std::size_t index_in_vector = std::distance(tlist.begin(),tit); 
//    const Triangle3 & t = *tit;
//    // dot product
//    double s = d * t.supporting_plane().orthogonal_vector();
//    if(s>0)
//    {
//      count += 1;
//    }else
//    {
//      count += -1;
//    }
//  }
//  w = count;
//
//}


void igl::winding_number_ray(
  const igl::embree::EmbreeIntersector & ei,
  const Eigen::MatrixXd & N,
  const Eigen::Vector3d & p, 
  const Eigen::MatrixXd & D,
  double & w,
  int & num_rays_shot)
{
  using namespace Eigen;
  // initialize count of total number of rays shot
  num_rays_shot = 0;
  // initialize winding number average
  w = 0;
  for(int s = 0;s<D.rows();s++)
  {
    Vector3d dir = D.row(s);
    int num_rays_shot_s = 0;
    double w_s = 0;
    winding_number_ray(ei,N,p,dir,w_s,num_rays_shot_s);
    w += w_s;
    num_rays_shot+=num_rays_shot_s;
  }
  // Compute average
  w = w/double(D.rows());
}

#ifdef WITH_CGAL
void winding_number_ray(
  const TriTree & tree,
  const TriangleList & tlist,
  const Eigen::Vector3d & p, 
  const Eigen::MatrixXd & D,
  double & w,
  int & num_rays_shot)
{
  using namespace Eigen;
  // initialize count of total number of rays shot
  num_rays_shot = 0;
  // initialize winding number average
  w = 0;
  for(int s = 0;s<D.rows();s++)
  {
    Vector3d dir = D.row(s);
    int num_rays_shot_s = 0;
    double w_s = 0;
    winding_number_ray(tree,tlist,p,dir,w_s,num_rays_shot_s);
    w += w_s;
    num_rays_shot+=num_rays_shot_s;
  }
  // Compute average
  w = w/double(D.rows());
}
#endif

#endif
