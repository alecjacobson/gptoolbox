#include "mex.h"
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/isolines.h>
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
  mexErrMsgTxt(nrhs>=4,"nrhs should be >= 4");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,F);
  Eigen::VectorXd S;
  parse_rhs_double(prhs+2,S);
  Eigen::VectorXd vals;
  parse_rhs_double(prhs+3,vals);

  Eigen::MatrixXd iV;
  Eigen::MatrixXi iE;
  Eigen::VectorXi I;

  igl::isolines(V,F,S,vals,iV,iE,I);

  switch(nlhs)
  {
    case 3:
      prepare_lhs_index(I,plhs+2);
    case 2:
      prepare_lhs_index(iE,plhs+1);
    case 1:
      prepare_lhs_double(iV,plhs+0);
    default:break;
  }
}


