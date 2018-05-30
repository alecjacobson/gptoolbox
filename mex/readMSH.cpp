#include "mex.h"

#include <igl/readMSH.h>
#include <igl/matlab/prepare_lhs.h>
#include <Eigen/Core>
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#  include <wordexp.h>
#endif


void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[])
{
  using namespace Eigen;
  /* Check for proper number of arguments */

  if (nrhs != 1) 
  {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
        "readMSH requires 1 input arguments, the path of the file to open");
  }

  // Read the file path
  char* file_path = mxArrayToString(prhs[0]);
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
  wordexp_t exp_result;
  wordexp(file_path, &exp_result, 0);
  file_path = exp_result.we_wordv[0];
#endif

  MatrixXd V;
  MatrixXi F;

  // Read the mesh
  if(!igl::readMSH(file_path,V,F))
  {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:fileio", "igl::readMSH failed.");
  }
  // Return the matrices to matlab
  switch(nlhs)
  {
    case 2:
      igl::matlab::prepare_lhs_index(F,plhs+1);
    case 1:
      igl::matlab::prepare_lhs_double(V,plhs);
    default: break;
  }

  return;
}


