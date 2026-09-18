#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/isolines_intrinsic.h>
#include <Eigen/Core>

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

  Eigen::MatrixXi F;
  mexErrMsgTxt(nrhs>=3,"nrhs should be >= 3");
  parse_rhs_index(prhs+0,F);
  Eigen::VectorXd S;
  parse_rhs_double(prhs+1,S);
  Eigen::VectorXd vals;
  parse_rhs_double(prhs+2,vals);

  Eigen::MatrixXd iB;
  Eigen::VectorXi iFI;
  Eigen::MatrixXi iE;
  Eigen::VectorXi I;

  igl::isolines_intrinsic(F,S,vals,iB,iFI,iE,I);

  switch(nlhs)
  {
    case 4:
      prepare_lhs_index(I,plhs+3);
    case 3:
      prepare_lhs_index(iE,plhs+2);
    case 2:
      prepare_lhs_index(iFI,plhs+1);
    case 1:
      prepare_lhs_double(iB,plhs+0);
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}
