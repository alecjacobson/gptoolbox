#include "mex.h"

#include <igl/read_triangle_mesh.h>
#include <igl/matlab/prepare_lhs.h>
#include <Eigen/Core>


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
        "readOBJ requires 1 input arguments, the path of the file to open");
  }

  // Read the file path
  char* file_path = mxArrayToString(prhs[0]);

  MatrixXd V;
  MatrixXi F;

  // Read the mesh
  if(!igl::read_triangle_mesh(file_path,V,F))
  {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:fileio", "igl::read_triangle_mesh failed.");
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

