#ifdef MEX
// mex -v -largeArrayDims -DMEX ...
//   -I/usr/local/igl/libigl/include ...
//   -I/opt/local/include/eigen3 ...
//   -I/opt/local/include ...
//   outer_hull.cpp

#include <igl/cgal/outer_hull.h>

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/C_STR.h>

#include <mex.h>
#include <Eigen/Dense>
#include <iostream>

#include <cstring>

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F)
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  mexErrMsgTxt(nrhs >= 2, "The number of input arguments must be >=2.");

  const auto & parse_mesh = [](
    const mxArray *prhs[], 
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
  {
    const int dim = mxGetN(prhs[0]);
    mexErrMsgTxt(dim == 3,
      "Mesh vertex list must be #V by 3 list of vertex positions");
    mexErrMsgTxt(dim == mxGetN(prhs[1]),
      "Mesh \"face\" simplex size must equal dimension");
    parse_rhs_double(prhs,V);
    parse_rhs_index(prhs+1,F);
  };
  parse_mesh(prhs,V,F);
}

void mexFunction(
  int nlhs, mxArray *plhs[], 
  int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;
  using namespace igl::cgal;

  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = cout.rdbuf(&mout);
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);

  MatrixXd V;
  MatrixXi F,G;
  VectorXi J,flip;
  parse_rhs(nrhs,prhs,V,F);
  outer_hull(V,F,G,J,flip);

  switch(nlhs)
  {
    default:
    {
      mexErrMsgTxt(false,"Too many output parameters.");
    }
    case 3:
    {
      prepare_lhs_logical(flip,plhs+2);
      // Fall through
    }
    case 2:
    {
      prepare_lhs_index(J,plhs+1);
      // Fall through
    }
    case 1:
    {
      prepare_lhs_index(G,plhs+0);
      // Fall through
    }
    case 0: break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

#endif


