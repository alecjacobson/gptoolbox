#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[])
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  mexErrMsgTxt(nrhs>=2,"nrhs should be == 2");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,F);

  switch(nlhs)
  {
    case 2:
      prepare_lhs_index(F,plhs+1);
    case 1:
      prepare_lhs_double(V,plhs+0);
    default:break;
  }
}

