#include <mex.h>
#include "segment_list.h"
#include "box_up.h"
#include <Eigen/Core>
#include <iostream>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/PI.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <igl/barycenter.h>
#include <algorithm>
namespace igl
{
// FORM_FACTOR
//
// [F,V,BC,N] = form_factor(P,E)
//
// Inputs:
//   P  #P by 2 list of 2D vertices
//   E  #E by 2 list of edge indices into P
// Outputs:
//   F  #E by #E matrix of integrated form factor values
//   V  #E by #E matrix of visibility values
//   BC  #E by 2 list of edge barycenters
//   N  #E by 2 list of length-weighted normals
template <
  typename DerivedP,
  typename DerivedE,
  typename DerivedF,
  typename DerivedV,
  typename DerivedBC,
  typename DerivedN>
inline void form_factor(
  const Eigen::MatrixBase<DerivedP> & P,
  const Eigen::MatrixBase<DerivedE> & E,
  Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedBC> & BC,
  Eigen::PlainObjectBase<DerivedN> & N)
{
  assert(P.cols() == 2 && "P should be 2D");
  assert(E.cols() == 2 && "E should be segments");
  const int m = E.rows();
  F.resize(m,m);
  V.resize(m,m);
  igl::barycenter(P,E,BC);
  DerivedN U;
  N.resize(m,2);
  for(int e = 0;e<m;e++)
  {
    N(e,0) = + (P(E(e,1),1)-P(E(e,0),1));
    N(e,1) = - (P(E(e,1),0)-P(E(e,0),0));
  }
  // unit normals
  U = N;
  U.rowwise().normalize();
  
  typedef int Index;
  typedef CGAL::Epick Kernel;
  // 2D Primitives
  typedef CGAL::Point_2<Kernel>    Point_2;
  typedef CGAL::Segment_2<Kernel>  Segment_2; 
  typedef std::vector<Segment_2> Segments;
  typedef typename std::vector<CGAL::Segment_2<Kernel> > Segments;
  typedef typename Segments::iterator SegmentsIterator;
  typedef typename Segments::const_iterator SegmentsConstIterator;
  typedef 
    CGAL::Box_intersection_d::Box_with_handle_d<double,2,SegmentsIterator> 
    Box;
  typedef 
    std::map<Index,std::vector<std::pair<Index,CGAL::Object> > > 
    OffendingMap;
  typedef std::map<std::pair<Index,Index>,std::vector<Index> >  EdgeMap;
  typedef std::pair<Index,Index> EMK;

  Segments scene_segments;
  std::vector<Box> scene_boxes;
  igl::copyleft::cgal::segment_list(P,E,scene_segments);

  igl::copyleft::cgal::box_up(scene_segments,scene_boxes);

  // rays as segments
  Segments ray_segments;
  ray_segments.reserve(m*(m-1)/2);
  std::vector<std::pair<int,int> > ray_subs;
  for(int i = 0;i<m;i++)
  {
    Point_2 pi(BC(i,0),BC(i,1));
    for(int j = i+1;j<m;j++)
    {
      Point_2 pj(BC(j,0),BC(j,1));
      ray_segments.emplace_back(pj,pi);
      ray_subs.push_back({i,j});
    }
  }
  assert(ray_segments.size() == (m*(m-1)/2));
  std::vector<Box> ray_boxes;
  igl::copyleft::cgal::box_up(ray_segments,ray_boxes);
  V.setConstant(m,m,true);
  const auto cb = [&](const Box &scene_box, const Box &ray_box) -> void
  {
    using namespace std;
    // index in F and T
    int scene_ind = scene_box.handle()-scene_segments.begin();
    int ray_ind = ray_box.handle()-ray_segments.begin();
    int source_ind = ray_subs[ray_ind].first;
    int dest_ind = ray_subs[ray_ind].second;
    if(source_ind == scene_ind || dest_ind == scene_ind)
    {
      return;
    }
    const Segment_2 & scene_seg = *scene_box.handle();
    const Segment_2 & ray_seg = *ray_box.handle();
    if(CGAL::do_intersect(scene_seg,ray_seg))
    {
      V(source_ind,dest_ind) = false;
      // also fill in transpose
      V(dest_ind,source_ind) = V(source_ind,dest_ind);
    }
  };
  CGAL::box_intersection_d(
    scene_boxes.begin(), scene_boxes.end(),
    ray_boxes.begin(), ray_boxes.end(),
    cb);
  typedef Eigen::Matrix<typename DerivedF::Scalar,1,2> RowVector2S;
  typedef typename DerivedF::Scalar Scalar;
  F.setConstant(m,m,0);
  for(int i = 0;i<m;i++)
  {
    const RowVector2S bci = BC.row(i);
    for(int j = i+1;j<m;j++)
    {
      const RowVector2S bcj = BC.row(j);
      const RowVector2S rij = bcj-bci;
      const Scalar rij_ui = +rij.dot(U.row(i));
      const Scalar rij_uj = -rij.dot(U.row(j));
      if( rij_ui > 0 && rij_uj > 0 && V(i,j))
      {
        const RowVector2S urij = rij.normalized();
        F(i,j) = (rij_ui*rij_uj)/(igl::PI*std::pow(rij.squaredNorm(),1.5));
        // also fill in transpose
        F(j,i) = F(i,j);
      }
    }
  }
}

}

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[]
)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  MatrixXd P;
  MatrixXi E;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs,P);
  parse_rhs_index(prhs+1,E);
  mexErrMsgTxt(P.cols()==2,"P must be #P by 2");
  mexErrMsgTxt(E.cols()==2,"E must be #E by 2");

  MatrixXd FF,V,BC,N;
  form_factor(P,E,FF,V,BC,N);

  switch(nlhs)
  {
    case 4:
      prepare_lhs_double(N,plhs+3);
    case 3:
      prepare_lhs_double(BC,plhs+2);
    case 2:
      prepare_lhs_double(V,plhs+1);
    case 1:
      prepare_lhs_double(FF,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}
