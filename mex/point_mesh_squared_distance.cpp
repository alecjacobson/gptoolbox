// 
//     mex -v -DMEX point_mesh_squared_distance.cpp ...
//       -Isrc -I/opt/local/include/ -I/opt/local/include/eigen3 ...
//       -I/usr/local/igl/libigl/include 
// 
#include <mex.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/point_mesh_squared_distance.h>

#ifdef MEX
#  include "mex.h"
#endif

#include <iostream>
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;

  MatrixXd P,V,C;
  VectorXi I;
  VectorXd sqrD;
  MatrixXi F;
  if(nrhs < 3)
  {
    mexErrMsgTxt("nrhs < 3");
  }
  parse_rhs_double(prhs,P);
  parse_rhs_double(prhs+1,V);
  parse_rhs_index(prhs+2,F);
  mexErrMsgTxt(P.cols()==3,"P must be #P by 3");
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3 || F.cols()==2 || F.cols()==1,"F must be #F by (3|2|1)");

  point_mesh_squared_distance(P,V,F,sqrD,I,C);
  // Prepare left-hand side
  switch(nlhs)
  {
    case 3:
    {
      // Treat indices as reals
      plhs[2] = mxCreateDoubleMatrix(C.rows(),C.cols(), mxREAL);
      double * Cp = mxGetPr(plhs[2]);
      copy(&C.data()[0],&C.data()[0]+C.size(),Cp);
      // Fallthrough
    }
    case 2:
    {
      // Treat indices as reals
      plhs[1] = mxCreateDoubleMatrix(I.rows(),I.cols(), mxREAL);
      double * Ip = mxGetPr(plhs[1]);
      VectorXd Id = (I.cast<double>().array()+1).matrix();
      copy(&Id.data()[0],&Id.data()[0]+Id.size(),Ip);
      // Fallthrough
    }
    case 1:
    {
      plhs[0] = mxCreateDoubleMatrix(sqrD.rows(),sqrD.cols(), mxREAL);
      double * sqrDp = mxGetPr(plhs[0]);
      copy(&sqrD.data()[0],&sqrD.data()[0]+sqrD.size(),sqrDp);
      break;
    }
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}
