#include "segment_list.h"
#include <Eigen/Core>
#include <iostream>
#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <algorithm>

// Inputs:
//   VA  #VA by 2 list of vertex positions
//   EA  #VA by 2 list of segment indices into VA
//   VB  #VB by 2 list of vertex positions
//   EB  #VB by 2 list of segment indices into VB
// Outputs:
//   IF  #IF by 2 list of pairs of indices into EA and EB
//   T  #IF by 2 list of paramteric distances along corresponding edge in EA and
//     EB. [-1 -1] indicates co-linear intersection.
void segment_segment_intersection(
  const Eigen::MatrixXd & VA,
  const Eigen::MatrixXi & EA,
  const Eigen::MatrixXd & VB,
  const Eigen::MatrixXi & EB,
  Eigen::MatrixXi & IF,
  Eigen::MatrixXd & T)
{

  typedef CGAL::Epeck Kernel;
  // 2D Primitives
  typedef CGAL::Point_2<Kernel>    Point_2;
  typedef CGAL::Segment_2<Kernel>  Segment_2; 
  std::vector<Segment_2> SA,SB;
  segment_list(VA,EA,SA);
  segment_list(VB,EB,SB);
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
  MatrixXd VA,VB;
  MatrixXi EA,EB;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  mexErrMsgTxt(nrhs==4,"nrhs should be == 4");
  parse_rhs_double(prhs,VA);
  parse_rhs_index(prhs+1,EA);
  mexErrMsgTxt(VA.cols()==2,"VA must be #VA by 2");
  mexErrMsgTxt(EA.cols()==2,"EA must be #EA by 2");
  parse_rhs_double(prhs+2,VB);
  parse_rhs_index(prhs+3,EB);
  mexErrMsgTxt(VB.cols()==2,"VB must be #VB by 2");
  mexErrMsgTxt(EB.cols()==2,"EB must be #EB by 2");

  switch(nlhs)
  {
    case 4:
    case 3:
    case 2:
    case 1:
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}
