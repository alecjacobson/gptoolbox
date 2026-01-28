#include <Eigen/Core>
#include <mex.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/C_STR.h>
#include <igl/matlab/requires_arg.h>
// i will be incremented on success
bool parse_spline_eytzinger_aabb(
  int & i,
  const int nrhs,
  const mxArray * prhs[],
  const char * name,
  Eigen::MatrixXd & B1,
  Eigen::MatrixXd & B2,
  Eigen::VectorXi & leaf)
{
  using namespace igl;
  using namespace igl::matlab;
  mexErrMsgTxt(mxIsCell(prhs[i+1]),
    C_STR("Parameter '"<<name<<"' requires Cell argument"));
  requires_arg(i,nrhs,name);
  const mxArray * cell = prhs[++i];
  mexErrMsgTxt(mxGetNumberOfElements(cell)==3,
    C_STR("Parameter '"<<name<<"' requires Cell argument with 3 elements"));
  // B1
  // ---
  {
    const mxArray * B1_mx = mxGetCell(cell,0);
    parse_rhs_double(&B1_mx,B1);
  }
  // B2
  // ---
  {
    const mxArray * B2_mx = mxGetCell(cell,1);
    parse_rhs_double(&B2_mx,B2);
  }
  // leaf
  // ---
  {
    const mxArray * leaf_mx = mxGetCell(cell,2);
    parse_rhs_index(&leaf_mx,leaf);
    mexErrMsgTxt(leaf.cols()==1, C_STR("leaf should have 1 column."));
  }
  return true;
}
