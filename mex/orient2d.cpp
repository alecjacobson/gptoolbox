#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/predicates/orient2d.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <Eigen/Core>

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[]
  )
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  mexErrMsgTxt(nrhs>=1,"nrhs should be == 1");
  Eigen::MatrixXd A,B,C,D;
  Eigen::VectorXi R;
  parse_rhs_double(prhs+0,A);
  parse_rhs_double(prhs+1,B);
  parse_rhs_double(prhs+2,C);
  igl::predicates::orient2d(A,B,C,R);

  switch(nlhs)
  {
    case 1:
      prepare_lhs_double(R,plhs+0);
    default:break;
  }
  std::cout.rdbuf(outbuf);
}


