#include <igl/simplify_polyhedron.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

#ifdef MEX
#  include <mex.h>
#  include <igl/C_STR.h>
#  include <igl/matlab/mexErrMsgTxt.h>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#endif

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>

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
  MatrixXd V,W;
  MatrixXi F,G;
  VectorXi J;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs,V);
  parse_rhs_index(prhs+1,F);
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
  simplify_polyhedron(V,F,W,G,J);

  switch(nlhs)
  {
    case 3:
      prepare_lhs_index(J,plhs+2);
    case 2:
      prepare_lhs_index(G,plhs+1);
    case 1:
      prepare_lhs_double(W,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}
