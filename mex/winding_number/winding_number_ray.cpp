#include "winding_number_ray.h"

#include <igl/embree/EmbreeIntersector.h>
#include <igl/embree/Hit.h>

#include <iterator> 
#include <iostream>

void igl::winding_number_ray_parity(
  const igl::EmbreeIntersector & ei,
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
  const igl::EmbreeIntersector & ei,
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
  const igl::EmbreeIntersector & ei,
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
